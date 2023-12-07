import os
import glob
import ntpath
import shutil
import string
import subprocess
from collections import Counter
import pandas as pd
import scipy.spatial
from Bio import SeqIO
from sklearn.cluster import DBSCAN
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
from networkx.algorithms.community import greedy_modularity_communities
from tqdm.auto import tqdm

# Dictionary mapping anticodon codes to values
ANTICODON_MAPPING = dict(GLYACC='gly4', HISATG='his2', CYSACA='cys2', ARGACG='arg1', ASPATC='asp2',
                         TYRATA='tyr2', PROAGG='pro4', ARGCCT='arg2', ALAAGC='ala3', SERAGA='ser6',
                         ASNATT='asn2', GLNCTG='gln2', GLUCTC='glu2', ARGCCG='arg4', THRAGT='thr4',
                         TRPCCA='trp1', GLYCCC='gly1', ILETAT='ile2', THRGGT='thr1', SERCGA='ser3',
                         ALACGC='ala4', PROCGG='pro1', PROGGG='pro2', LEUTAG='leu3', SERGGA='ser2',
                         LEUTAA='leu2', ALAGGC='ala1', VALTAC='val1', THRCGT='thr2', TYRGTA='tyr1',
                         ASPGTC='asp1', HISGTG='his1', GLUTTC='glu1', ASNGTT='asn1', GLNTTG='gln1',
                         LEUAAG='leu6', PHEAAA='phe2', VALAAC='val4', LYSTTT='lys1', METCAT='met1',
                         ARGTCT='arg3', ILEGAT='ile1', LYSCTT='lys2', ALATGC='ala2', SERTGA='ser1',
                         PROTGG='pro3', LEUGAG='leu1', ARGTCG='arg5', VALGAC='val2', GLYTCC='gly3',
                         PHEGAA='phe1', SECTCA='sec1', CYSGCA='cys1', GLYGCC='gly2', ARGGCG='arg6',
                         ILEAAT='ile3', VALCAC='val3', LEUCAA='leu5', LEUCAG='leu4', THRTGT='thr3',
                         SERGCT='ser4', SERACT='ser5', SUPCTA='Sup1', SUPTTA='Sup2', SUPTCA='Sup3',
                         ILE2CAT='ile3', FMETCAT='fmet1', ssrA='ssrA')


def apply_anticodon_mapping(id_info):
    if "Name=tmRNA" in id_info:
        return "ssrA"
    else:
        anticodon = id_info.split("Name=")[1].replace("-", "").upper()
        return ANTICODON_MAPPING.get(anticodon, '')


def extract_accession_dict(list_files):
    accession_dict = {}
    for in_file in list_files:
        with open(in_file, "r") as strain_file:
            for record in SeqIO.parse(in_file, "fasta"):
                accession_dict[record.id] = ntpath.basename(in_file).replace(f".fasta", "")

    return accession_dict


def check_core(no_strains, lcb_group):
    codes = [int(i.split("\t")[0]) for i in lcb_group[2:]]
    for num in range(1, no_strains + 1):
        if num not in codes:
            return False
    return True


def create_dict_row(block, dct_codes, dct_files):
    new_rows = []
    ident = block[0].split()[1][1:]
    for rline in block[2:]:
        lcb = rline.strip().split("\t")
        if int(lcb[-3]) < int(lcb[-2]):
            lcb_start = int(lcb[-3])
            lcb_end = int(lcb[-2])
        else:
            lcb_start = int(lcb[-2])
            lcb_end = int(lcb[-3])
        new_row = {
            "strain": dct_files[dct_codes[lcb[0]]],
            "method": "SibeliaZ-core",
            "type": "LCB-Core",
            "start": lcb_start,
            "end": lcb_end,
            "na1": ".",
            "sense": lcb[-4],
            "na2": ".",
            "ID": f"Name=Block_#{ident}"
        }
        new_rows.append(new_row)
    return new_rows


def create_lcbs_df(blocks_file, list_files):
    df = pd.DataFrame(columns=["strain", "method", "type", "start", "end", "na1", "sense", "na2", "ID"])
    dct_files = extract_accession_dict(list_files)
    with open(blocks_file, "r") as lcbs_file:
        fields = "".join(lcbs_file.readlines()).split(f"{80 * '-'}\n")
        seqs = [i.strip().split("\n") for i in fields if not "Block #" in i][0][1:]
        groups = [i.strip().split("\n") for i in fields if "Block #" in i]
        conserved_lcb = [i for i in groups]
    dct_codes = {}
    for i in seqs:
        line = i.strip().split("\t")
        dct_codes[line[0]] = line[-1]
    rows_lcb = [create_dict_row(block, dct_codes, dct_files) for block in conserved_lcb]
    new_rows = [row for group in rows_lcb for row in group]
    df_new = pd.DataFrame(new_rows, columns=["strain", "method", "type", "start", "end", "na1", "sense", "na2", "ID"])
    df = pd.concat([df, df_new])
    for group_name, df_group in df.groupby("strain"):
        df_group['ID'] = df_group['ID'] + '_' + df_group.groupby('ID').cumcount().add(1).astype(str)
        df.update(df_group)
    return df


def create_tmdnas_df(list_files):
    df = pd.DataFrame(columns=["strain","method","type","start","end","na1","sense","na2","ID"])
    for infile in list_files:
        df2 = pd.read_csv(infile, sep="\t", names=["strain","method","type","start","end","na1","sense","na2","ID"])
        df = pd.concat([df,df2], ignore_index=True)
    return df


def delete_non_core(lcbdf):
    groups = lcbdf.groupby("ID", group_keys=False)
    filtered_groups = groups.filter(lambda x: len(x) >= 0.8 * lcbdf['strain'].nunique())
    new_df = filtered_groups.groupby("ID", group_keys=False).apply(lambda x: pd.concat([x]))
    new_df.rename(columns={'ID': 'IDs'}, inplace=True)
    return new_df


def run_sibeliaz(list_files, output_folder, threads):
    if not os.path.exists(f'{output_folder}/all_chr_sibelia/blocks_coords.gff'):
        cmd_sibeliaz = f"sibeliaz -k 11 -n -t {threads} -o {output_folder}/all_chr_sibelia/ {' '.join(list_files)} ; maf2synteny -b 1000 -o {output_folder}/all_chr_sibelia/ {output_folder}/all_chr_sibelia/blocks_coords.gff"
        p = subprocess.run(cmd_sibeliaz, shell=True)


def seqlen_dct(list_files, file_ext):
    str_dct = {}
    for in_file in list_files:
        with open(in_file, "r") as strain_file:
            for record in SeqIO.parse(in_file, "fasta"):
                str_dct[ntpath.basename(in_file).replace(f".{file_ext}", "")] = len(record.seq)
    return str_dct


def overlap_at_start(r1, r2):
    return r1[0] == r2[0] or (r1[0] < r2[0] and r1[1] > r2[0])


def overlap_at_end(r1, r2):
    return r1[1] == r2[1] or (r1[0] < r2[1] and r1[1] > r2[1])


def check_if_overlap(range1, range2):
    return not (range1[1] < range2[0] or range2[1] < range1[0])


def find_nearest_range(point, ranges):
    nearest_range = min(ranges, key=lambda r: abs(point - r[0]))
    return nearest_range


def up_lcb_finder(row, dct_len, lcbdf):
    len_seq = dct_len[row.strain]
    strain_lcb = lcbdf[lcbdf['strain'] == row.strain]
    strain_lcb = strain_lcb.sort_values(by=["start"])
    lcb_tuples_ori = [(start, end) for start, end in zip(strain_lcb['start'], strain_lcb['end'])]
    lcb_tuples = [(start, end) for start, end in zip(strain_lcb['start'], strain_lcb['end'])]
    lcb_tuples.append((lcb_tuples[-1][0] - len_seq, lcb_tuples[-1][1] - len_seq))
    lcb_tuples.append((lcb_tuples[0][0] + len_seq, lcb_tuples[0][1] + len_seq))
    non_overlap = [lcb for lcb in lcb_tuples if not check_if_overlap(lcb, (row.start, row.end))]
    lcb_overlap = [lcb for lcb in lcb_tuples if overlap_at_start(lcb, (row.start, row.end))]
    nearest_ur = find_nearest_range(row.start, [ranges for ranges in non_overlap if ranges[0] < row.start])
    if len(lcb_overlap) != 0:
        prev_start = lcb_overlap[0][0]
    else:
        prev_start = nearest_ur[0]
    if prev_start < 0:
        lcb_start = lcb_tuples_ori[-1][0]
    elif prev_start > len_seq:
        lcb_start = lcb_tuples_ori[0][0]
    else:
        lcb_start = prev_start
    lcb_ur = strain_lcb.loc[strain_lcb["start"] == lcb_start]
    return((lcb_ur.iloc[0]["IDs"], lcb_ur.iloc[0]["start"], lcb_ur.iloc[0]["end"], lcb_ur.iloc[0]["sense"]))

def down_lcb_finder(row, dct_len, lcbdf):
    len_seq = dct_len[row.strain]
    strain_lcb = lcbdf[lcbdf['strain'] == row.strain]
    strain_lcb = strain_lcb.sort_values(by=["start"])
    lcb_tuples_ori = [(start, end) for start, end in zip(strain_lcb['start'], strain_lcb['end'])]
    lcb_tuples = [(start, end) for start, end in zip(strain_lcb['start'], strain_lcb['end'])]
    lcb_tuples.append((lcb_tuples[-1][0] - len_seq, lcb_tuples[-1][1] - len_seq))
    lcb_tuples.append((lcb_tuples[0][0] + len_seq, lcb_tuples[0][1] + len_seq))
    non_overlap = [lcb for lcb in lcb_tuples if not check_if_overlap(lcb, (row.start, row.end))]
    lcb_overlap = [lcb for lcb in lcb_tuples if overlap_at_end(lcb, (row.start, row.end))]
    nearest_ur = find_nearest_range(row.start, [ranges for ranges in non_overlap if ranges[0] > row.end])
    if len(lcb_overlap) != 0:
        prev_start = lcb_overlap[0][0]
    else:
        prev_start = nearest_ur[0]
    if prev_start < 0:
        lcb_start = lcb_tuples_ori[-1][0]
    elif prev_start > len_seq:
        lcb_start = lcb_tuples_ori[0][0]
    else:
        lcb_start = prev_start
    lcb_ur = strain_lcb.loc[strain_lcb["start"] == lcb_start]
    return((lcb_ur.iloc[0]["IDs"], lcb_ur.iloc[0]["start"], lcb_ur.iloc[0]["end"], lcb_ur.iloc[0]["sense"]))


def UR_block_name(row):
    if row["sense"] == "+":
        return row["LCB_UP"][0].replace("Name=Block_", "")
    else:
        return row["LCB_DOWN"][0].replace("Name=Block_", "")


def DR_block_name(row):
    if row["sense"] == "+":
        return row["LCB_DOWN"][0].replace("Name=Block_", "")
    else:
        return row["LCB_UP"][0].replace("Name=Block_", "")


def tdnas_neigh_finder(row, dct_len, tmrnadf):
    len_seq = dct_len[row.strain]
    if row.sense == "+":
        range_neigh = (row.start - 1000, row.start - 1)
    else:
        range_neigh = (row.end + 1, row.end + 1000)
    strain_tmrnas = tmrnadf[tmrnadf['strain'] == row.strain]
    strain_tmrnas = strain_tmrnas.sort_values(by=["start"])
    tmrnas_tuples = [(start, end) for start, end in zip(strain_tmrnas['start'], strain_tmrnas['end'])]
    tmrnas_neigh = [tmdna for tmdna in tmrnas_tuples if check_if_overlap(tmdna, range_neigh)]
    anticodon = []
    for tmrna in tmrnas_neigh:
        row = tmrnadf.loc[tmrnadf["start"] == tmrna[0]]
        anticodon.append(row.iloc[0]["ANT"])
    return "-".join(anticodon)


def check_UR_sense(row):
    if row.sense == "+":
        LCB_UR = row.LCB_UP
    else:
        LCB_UR = row.LCB_DOWN
    if LCB_UR[3] == row.sense:
        return LCB_UR[0] + "_YES"
    else:
        return LCB_UR[0] + "_NO"


def create_dict_ctxs(df):
    dict_ctx = {}
    for group_name, df_group in df.groupby("ANT"):
        if len(df_group) > 1:
            df_clust = df_group[["UR_name", "tdnas_neigh"]]
            df_clust_binary = pd.get_dummies(df_clust)
            jaccard = scipy.spatial.distance.cdist(df_clust_binary, df_clust_binary, metric='jaccard')
            rows_distance = pd.DataFrame(jaccard, columns=df_clust_binary.index.values,
                                         index=df_clust_binary.index.values)
            # Instantiate the DBSCAN object with the desired hyperparameters
            G = nx.Graph()
            G.add_nodes_from(df.columns)
            for i, row in rows_distance.iterrows():
                for j, value in row.items():
                    if i != j and value != 1.0:
                        G.add_edge(i, j, weight=value)
            labels = list(louvain_communities(G))
            # Perform DBSCAN clustering on the distance matrix
        else:
            labels = [0]

        # assigns new clusters for ourliers
        for i in range(len(labels)):
            if labels[i] < 0:
                labels[i] = max(labels) + abs(labels[i])

        # alphabet list
        alphabet_list = list(string.ascii_uppercase)
        double_letter_list = [letter1 + letter2 for letter1 in alphabet_list for letter2 in alphabet_list]
        alphabet = alphabet_list + double_letter_list

        # create a list of clusters and sort based in frequency (most common first)
        clusters = list(dict.fromkeys([item for items, c in Counter(labels).most_common() for item in [items] * c]))

        code_dict = {}
        for index, value in enumerate(clusters):
            code_dict[value] = alphabet[index]

        letter_dict = {}
        x = 0
        for index, row in df_group.iterrows():
            letter_dict[row.name] = code_dict[labels[x]]
            x = x + 1

        dict_ctx.update(letter_dict)
    return dict_ctx


def apply_nom(row, dict_ctx, nom_ext):
    tdna_name = f"{nom_ext}_{row.ANT}{str(dict_ctx[row.name])}"
    return tdna_name


def DR_block_name(row):
    if row["sense"] == "+":
        return row["LCB_DOWN"][0].replace("Name=Block_","")
    else:
        return row["LCB_UP"][0].replace("Name=Block_","")


def tdna_clusterization(input_folder, output_folder, threads, prefix):
    # Step 1: Run SibeliaZ to identify LCBs and create synteny blocks
    list_files = glob.glob(f'{output_folder}/*/*.fasta')
    run_sibeliaz(list_files, output_folder, threads)

    # Step 2: Create DataFrames for LCBs and tDNAs
    blocks_file = f'{output_folder}/all_chr_sibelia/1000/blocks_coords.txt'
    lcbdf = create_lcbs_df(blocks_file, list_files)

    tdnas_files = glob.glob(f'{output_folder}/*/*.tmdnas')
    tdnadf = create_tmdnas_df(tdnas_files)

    # Step 3: Delete non-core LCBs
    core_lcbdf = delete_non_core(lcbdf)

    # Step 4: Find nearest upstream and downstream LCBs for each tDNA
    dct_len = seqlen_dct(list_files, "fasta")
    tdnadf["LCB_UP"] = tdnadf.apply(up_lcb_finder, args=(dct_len, core_lcbdf), axis=1)
    tdnadf["LCB_DOWN"] = tdnadf.apply(down_lcb_finder, args=(dct_len, core_lcbdf), axis=1)

    # Step 5: Find anticodon for each tDNA
    tdnadf["ANT"] = tdnadf["ID"].apply(apply_anticodon_mapping)

    # Step 6: Find neighboring tDNAs for each tDNA
    tdnadf["tdnas_neigh"] = tdnadf.apply(tdnas_neigh_finder, args=(dct_len, tdnadf), axis=1)

    # Step 7: Create context clusters based on tDNA neighbors
    tdnadf["UR_name"] = tdnadf.apply(UR_block_name, axis=1)
    tdnadf["DR_name"] = tdnadf.apply(lambda row: DR_block_name(row), axis=1)
    dict_ctx = create_dict_ctxs(tdnadf)

    tdnadf["tdna_name"] = tdnadf.apply(lambda row: apply_nom(row, dict_ctx, prefix), axis=1)
    tdnadf[["strain", "start", "end", "sense", "LCB_UP", "LCB_DOWN","UR_name", "DR_name", "tdna_name"]].to_csv(
        f"{output_folder}/{prefix}_tdna_scheme.csv")
    lcbdf.to_csv(f"{output_folder}/{prefix}_LCBs.csv")


