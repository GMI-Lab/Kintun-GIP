import glob
from Bio import SeqIO
import os
import pandas as pd
from tqdm import tqdm
import ast
import sys

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
    # tRNA-Ala(ggc)
    if "transfer-messenger" in id_info:
        return "ssrA"
    else:
        anticodon = id_info.split("-")[1].replace("(", "").replace(")", "").upper()
        return ANTICODON_MAPPING.get(anticodon, '')

def set_markers(row):
    if row.Strand == 1:
        ctx_code = set(row.UP_neigh) | set([row.tdnas_neigh])
    else:
        ctx_code = set(row.DOWN_neigh) | set([row.tdnas_neigh])
    return ",".join(sorted([str(i) for i in list(ctx_code)]))

# Get unique values and their frequencies, then sort by frequency
def generate_label(n, anticodon, prefix):
    # Generate alphabetical labels (A, B, ..., Z, AA, AB, ..., ZZ, AAA, AAB, ...)
    result = ""
    while n >= 0:
        result = chr(n % 26 + ord('A')) + result
        n = n // 26 - 1
    return f"{prefix}_" + anticodon + result


def find_neighbors(tRNA_row, cds_df, gene_ids):
    # Filter CDS rows based on the same "Strand"
    neighbors_cds = cds_df[cds_df['File'] == tRNA_row['File']]

    # Filter CDS rows based on their position relative to tRNA
    up_neighbors = neighbors_cds[
        (neighbors_cds['End'] < tRNA_row['Start']) |
        ((neighbors_cds['Start'] > tRNA_row['End']) & (neighbors_cds['Start'] <= tRNA_row['Start']))
    ].tail(3)
    if len(up_neighbors) < 2:
        up_neighbors = neighbors_cds.head(3)

    down_neighbors = neighbors_cds[
        (neighbors_cds['Start'] > tRNA_row['End']) |
        ((neighbors_cds['End'] < tRNA_row['Start']) & (neighbors_cds['End'] >= tRNA_row['End']))
    ].head(3)
    if len(down_neighbors) < 2:
        down_neighbors = neighbors_cds.head(3)


    up_neigh = []
    for index, row in up_neighbors.iterrows():
        if tRNA_row.Strand == row.Strand:
            up_neigh.append(f"{gene_ids[row.ID]}_Yes")
        else:
            up_neigh.append(f"{gene_ids[row.ID]}_No")

    down_neigh = []
    for index, row in down_neighbors.iterrows():
        if tRNA_row.Strand == row.Strand:
            down_neigh.append(f"{gene_ids[row.ID]}_Yes")
        else:
            down_neigh.append(f"{gene_ids[row.ID]}_No")

    up_ids = up_neighbors[['ID']].values.tolist() if not up_neighbors.empty else []
    down_ids = down_neighbors[['ID']].values.tolist() if not down_neighbors.empty else []

    return (up_neigh, down_neigh, up_ids, down_ids)



def parse_genbank_files_to_dataframe(genbank_files):
    data = {
        'File': [],
        'Chr_size': [],
        'ID': [],
        'Type': [],
        'Product': [],
        'Start': [],
        'End': [],
        'Strand': []
    }

    for genbank_file in genbank_files:
        records = SeqIO.parse(genbank_file, "genbank")

        for record in records:
            for feature in record.features:
                if feature.type == 'tRNA' or feature.type == 'tmRNA':
                    data['File'].append(os.path.basename(genbank_file).replace(".gbff", ""))
                    data['Chr_size'].append(len(record.seq))
                    data['ID'].append(feature.qualifiers.get('locus_tag', [None])[0])
                    data['Type'].append(feature.type)
                    data['Product'].append(feature.qualifiers.get('product', [None])[0])
                    data['Start'].append(feature.location.start)
                    data['End'].append(feature.location.end)
                    data['Strand'].append(feature.location.strand)

    df = pd.DataFrame(data)
    return df


def rrna_parse_genbank_files_to_dataframe(genbank_files):
    data = {
        'File': [],
        'Chr_size': [],
        'ID': [],
        'Type': [],
        'Product': [],
        'Start': [],
        'End': [],
        'Strand': []
    }

    for genbank_file in tqdm(genbank_files, desc="Processing rDNAs in GenBank files"):
        records = SeqIO.parse(genbank_file, "genbank")

        for record in records:
            for feature in record.features:
                if feature.type == 'rRNA':
                    data['File'].append(os.path.basename(genbank_file).replace(".gbff", ""))
                    data['Chr_size'].append(len(record.seq))
                    data['ID'].append(feature.qualifiers.get('locus_tag', [None])[0])
                    data['Type'].append(feature.type)
                    data['Product'].append(feature.qualifiers.get('product', [None])[0])
                    data['Start'].append(feature.location.start)
                    data['End'].append(feature.location.end)
                    data['Strand'].append(feature.location.strand)

    df = pd.DataFrame(data)
    return df


def cds_parse_genbank_files_to_dataframe(genbank_files, core_tags):
    data = {
        'File': [],
        'ID': [],
        'Type': [],
        'Start': [],
        'End': [],
        'Strand': []
    }

    for genbank_file in tqdm(genbank_files, desc="Processing CDSs in GenBank files"):
        records = SeqIO.parse(genbank_file, "genbank")
        for record in records:
            for feature in [feature for feature in record.features if
                            [feature.qualifiers.get('locus_tag', [None])[0]] in core_tags[
                                os.path.basename(genbank_file).replace(".gbff", "")] and feature.type == "CDS"]:
                data['File'].append(os.path.basename(genbank_file).replace(".gbff", ""))
                data['ID'].append(feature.qualifiers.get('locus_tag', [None])[0])
                data['Type'].append(feature.type)
                data['Start'].append(feature.location.start)
                data['End'].append(feature.location.end)
                data['Strand'].append(feature.location.strand)

    df = pd.DataFrame(data)
    return df


def check_if_overlap(range1, range2):
    return not (range1[1] < range2[0] or range2[1] < range1[0])


def tdnas_neigh_finder(row, tmrnadf):
    ret_res = False
    if row.Strand == 1:
        range_neigh = (row.Start - 6000, row.Start - 1)
    else:
        ret_res = True
        range_neigh = (row.End + 1, row.End + 6000)
    strain_tmrnas = tmrnadf[tmrnadf['File'] == row.File]
    strain_tmrnas = strain_tmrnas.sort_values(by=["Start"])
    tmrnas_tuples = [(start, end) for start, end in zip(strain_tmrnas['Start'], strain_tmrnas['End'])]
    tmrnas_neigh = [tmdna for tmdna in tmrnas_tuples if check_if_overlap(tmdna, range_neigh)]
    anticodon = []
    for tmrna in tmrnas_neigh:
        row = tmrnadf.loc[tmrnadf["Start"] == tmrna[0]]
        anticodon.append(row.iloc[0]["ANT"])
    if ret_res:
        return "-".join(anticodon)
    else:
        return "-".join(anticodon[::-1])


def check_if_overlap(range1, range2):
    return not (range1[1] < range2[0] or range2[1] < range1[0])


def check_if_overlap_circular(range1, range2, chromosome_size):
    if check_if_overlap(range1, range2):
        return True
    # Check for circular overlap in both directions
    circular_range_forward = (range2[0] + chromosome_size, range2[1] + chromosome_size)
    circular_range_backward = (range2[0] - chromosome_size, range2[1] - chromosome_size)

    if check_if_overlap(circular_range_forward, range1) or check_if_overlap(circular_range_backward, range1):
        return True
    return False


def check_cds(tdnas_df, cds_df, core_df, cds_df_dict):
    left_dict = {}
    right_dict = {}
    tdnas_groups = tqdm(tdnas_df.groupby("tdna_name"), desc="Processing tDNAs", total=len(tdnas_df["tdna_name"].unique()))
    for group_name, group_df in tdnas_groups:
        LEFT_groups_dict = {}
        RIGHT_groups_dict = {}
        LEFT_groups = []
        RIGHT_groups = []
        for index, row in group_df.iterrows():
            for index,value in enumerate(row.UP_neigh):
                LEFT_groups.append(value.replace("_Yes", "").replace("_No", ""))
                LEFT_groups_dict[(value.replace("_Yes", "").replace("_No", ""), row.File)] = row.UP_neigh_ori[index][0]
            for index,value in enumerate(row.DOWN_neigh):
                RIGHT_groups.append(value.replace("_Yes", "").replace("_No", ""))
                RIGHT_groups_dict[(value.replace("_Yes", "").replace("_No", ""), row.File)] = row.DOWN_neigh_ori[index][0]
        LEFT_groups = list(set(LEFT_groups))
        RIGHT_groups = list(set(RIGHT_groups))
        for index, row in group_df.iterrows():
            LEFT_coords = []
            RIGHT_coords = []
            for k in LEFT_groups:
                try:
                    LEFT_coords.append(cds_df_dict[(LEFT_groups_dict[(k, row.File)], row.File)])
                except KeyError:
                    continue
            for l in RIGHT_groups:
                try:
                    RIGHT_coords.append(cds_df_dict[(RIGHT_groups_dict[(l, row.File)], row.File)])
                except KeyError:
                    continue
            ctx_range = (row.Start - 800000, row.End + 800000)
            try:
                left_dict[row.name] = (min([range1 for range1 in LEFT_coords if check_if_overlap_circular(range1, ctx_range, row.Chr_size)]))
            except ValueError:
               print(row)
               sys.exit()
            try:
                right_dict[row.name] = (max([range1 for range1 in RIGHT_coords if check_if_overlap_circular(range1, ctx_range, row.Chr_size)]))
            except ValueError:
                print(row)
                sys.exit()
    return right_dict, left_dict


def check_overlap(conserved_blocks, dict_clusters, tdna_name):
    ur_dict = {}
    dr_dict = {}
    for value, ranges in conserved_blocks.items():
        target_start = int(value.split("_")[2])
        target_end = int(value.split("_")[3])
        updated_ranges = []
        for frag_start, frag_end in ranges:
            frag_start = int(frag_start)
            frag_end = int(frag_end)
            if frag_end < target_start:
                updated_ranges.append((frag_start, frag_end))
            elif frag_start > target_end:
                updated_ranges.append((frag_start, frag_end))
            else:
                if frag_start < target_start and frag_end > target_end:
                    updated_ranges.append((frag_start, target_start - 1))
                    if frag_end - target_end + 1 > 500:
                        updated_ranges.append((target_end + 1, frag_end))
                elif frag_end >= target_start and frag_end <= target_end:
                    updated_ranges.append((frag_start, target_start - 1))
                elif frag_start >= target_start and frag_start <= target_end:
                    updated_ranges.append((target_end + 1, frag_end))
            updated_ranges = sorted(updated_ranges, key=lambda x: x[0])
        for frag_start, frag_end in updated_ranges:
            if frag_end < target_start:
                ur_dict[value] = (frag_start, frag_end)
            elif frag_start > target_end:
                dr_dict[value] = (frag_start, frag_end)
                break
    return ur_dict, dr_dict


def circular_range_length(pilr_range, chromosome_length):
    start = pilr_range[0]
    end = pilr_range[1]
    if start <= end:
        return min(end - start, chromosome_length - (end - start))
    else:
        return min(start - end, chromosome_length - (start - end))


def PILR_prediction(row, value):
    if value == "down":
        if row.Strand == 1:
            pilr_range = (row.End + 1, row.block_right[0] - 1)
        else:
            pilr_range = (row.block_left[1] + 1, row.Start -1)
    if value == "up":
        if row.Strand == -1:
            pilr_range = (row.End + 1, row.block_right[0] - 1)
        else:
            pilr_range = (row.block_left[1] + 1, row.Start - 1)
    if circular_range_length(pilr_range, row.Chr_size) >= 5000:
        return pilr_range
    else:
        return "None"


def get_range(row, dict_pilrs):
    return dict_pilrs[row.name]


def safe_literal_eval(cell):
    try:
        return ast.literal_eval(cell)
    except (SyntaxError, ValueError):
        return cell


def tdna_clusterization(input_folder, output_folder, prefix, format):
    # Parse a list of GenBank files and create a DataFrame
    genbank_files_list = glob.glob(f"{input_folder}/*.{format}")
    tdnas_df = parse_genbank_files_to_dataframe(genbank_files_list)
    rdnas_df = rrna_parse_genbank_files_to_dataframe(genbank_files_list)
    rdnas_df.to_csv(f"{output_folder}/{prefix}_rDNAs.csv")

    # Decorate pandas' apply functions with tqdm
    tqdm.pandas()

    # Accelerate the processing of pd.read_csv with tqdm
    panaroo_roary = pd.read_csv(f"{input_folder}/gene_presence_absence_roary.csv", header=0, low_memory=False).progress_apply(lambda x: x)
    panaroo_roary = panaroo_roary.loc[(panaroo_roary["No. isolates"] > 0.9 * len(genbank_files_list)) & (
                panaroo_roary["No. sequences"] > 0.9 * len(genbank_files_list))]

    panaroo_roary.to_csv(f"{output_folder}/{prefix}_coregenome_detail.csv")

    # Create a dictionary
    core_tags = {}

    # Use tqdm to iterate over columns, starting from the 15th column
    for column in tqdm(panaroo_roary.columns[14:], desc="Processing columns"):
        # Split values by semicolon if present, and convert to list
        values = panaroo_roary[column].apply(lambda x: x.split(';') if ';' in str(x) else [x]).tolist()
        # Store values in the core_tags dictionary
        core_tags[column] = values

    gene_ids = {}
    # Use tqdm to iterate over rows of the DataFrame
    for index, row in tqdm(panaroo_roary.iterrows(), total=len(panaroo_roary), desc="Processing rows"):
        key = row['Gene']
        # Iterate over columns starting from the 15th column
        for column, value in row.iloc[14:].items():
            # Split values by semicolon if present
            value_list = str(value).split(';')
            for split_value in value_list:
                gene_ids[split_value] = key

    cds_df = cds_parse_genbank_files_to_dataframe(genbank_files_list, core_tags)
    cds_df.to_csv(f"{output_folder}/{prefix}_core_CDSs.csv")

    tdnas_df["ANT"] = tdnas_df["Product"].apply(apply_anticodon_mapping)
    tdnas_df = tdnas_df[tdnas_df['Product'] != 'tRNA-Xxx']
    neighbors_series = tdnas_df.apply(lambda row: find_neighbors(row, cds_df, gene_ids), axis=1)
    tdnas_df['UP_neigh'] = neighbors_series.apply(lambda x: x[0])
    tdnas_df['DOWN_neigh'] = neighbors_series.apply(lambda x: x[1])
    tdnas_df['UP_neigh_ori'] = neighbors_series.apply(lambda x: x[2])
    tdnas_df['DOWN_neigh_ori'] = neighbors_series.apply(lambda x: x[3])
    tdnas_df["tdnas_neigh"] = tdnas_df.apply(tdnas_neigh_finder, args=(tdnas_df,), axis=1)

    tdnas_df = tdnas_df.map(safe_literal_eval)

    tdnas_df["up_code"] = tdnas_df.apply(set_markers, axis=1)

    # try to mapping 
    mapping = {}
    for group_name, group_df in tqdm(tdnas_df.groupby("ANT"), desc="Mapping groups"):
        # Sort values and count occurrences
        sorted_counts = group_df['up_code'].value_counts().sort_values(ascending=False)
        # Generate labels and update the mapping dictionary
        part_mapping = {value: generate_label(i, group_name, prefix) for i, value in enumerate(sorted_counts.index)}
        mapping.update(part_mapping)

    tdnas_df['tdna_name'] = tdnas_df['up_code'].map(mapping)

    # check left and right ctx

    cds_df_dict = {}
    for index, row in cds_df.iterrows():
        cds_df_dict[(row.ID,row.File)] = (row.Start,row.End)
    right_blocks, left_blocks = check_cds(tdnas_df, cds_df, panaroo_roary, cds_df_dict)
    tdnas_df["block_left"] = tdnas_df.apply(get_range, args=(left_blocks,), axis=1)
    tdnas_df["block_right"] = tdnas_df.apply(get_range, args=(right_blocks,), axis=1)

    # check for upstream and downstream pilrs
    tdnas_df["PILR_UR"] = tdnas_df.apply(PILR_prediction, args=("up",), axis=1)
    tdnas_df["PILR_DR"] = tdnas_df.apply(PILR_prediction, args=("down",), axis=1)

    # Calculate Prevalence
    value_counts = tdnas_df["tdna_name"].value_counts()
    n_strain = tdnas_df["File"].nunique()

    # Calculate prevalence (frequency) for each value
    prevalence = tdnas_df['tdna_name'].map(lambda x: value_counts[x] / n_strain)

    # Add the prevalence column to the dataframe
    tdnas_df['Prevalence'] = prevalence

    for index, row in tdnas_df.iterrows():
        if row["PILR_DR"] != "None":
            start = row['PILR_DR'][0]
            end = row['PILR_DR'][1]
            file = row["File"]
            other_tRNAs = tdnas_df[(tdnas_df['Start'] >= start) & (tdnas_df['End'] <= end) & (tdnas_df.index != index) & (tdnas_df.File == file)]
            # Check if other tRNA genes have a prevalence over 0.8
            high_prevalence_tRNAs = other_tRNAs[other_tRNAs['Prevalence'] > 0.8]
            # If high prevalence tRNA genes are found, update the range to zero
            if not high_prevalence_tRNAs.empty:
                tdnas_df.at[index, 'PILR_DR'] = "None"

    for index, row in tdnas_df.iterrows():
        if row["PILR_UR"] != "None":
            start = row['PILR_UR'][0]
            end = row['PILR_UR'][1]
            file = row["File"]
            other_tRNAs = tdnas_df[(tdnas_df['Start'] >= start) & (tdnas_df['End'] <= end) & (tdnas_df.index != index) & (tdnas_df.File == file)]
            # Check if other tRNA genes have a prevalence over 0.8
            high_prevalence_tRNAs = other_tRNAs[other_tRNAs['Prevalence'] > 0.8]
            # If high prevalence tRNA genes are found, update the range to zero
            if not high_prevalence_tRNAs.empty:
                tdnas_df.at[index, 'PILR_UR'] = "None"

    for index, row in tdnas_df.iterrows():
        if row["PILR_DR"] != "None":
            start = row['PILR_DR'][0]
            end = row['PILR_DR'][1]
            file = row["File"]
            other_rRNAs = rdnas_df[(rdnas_df['Start'] >= start) & (rdnas_df['End'] <= end)  & (rdnas_df.File == file)]
            if not other_rRNAs.empty:
                tdnas_df.at[index, 'PILR_DR'] = "None"

    for index, row in tdnas_df.iterrows():
        if row["PILR_UR"] != "None":
            start = row['PILR_UR'][0]
            end = row['PILR_UR'][1]
            file = row["File"]
            other_rRNAs = rdnas_df[(rdnas_df['Start'] >= start) & (rdnas_df['End'] <= end)  & (rdnas_df.File == file)]
            if not other_rRNAs.empty:
                tdnas_df.at[index, 'PILR_DR'] = "None"

    tdnas_df.to_csv(f"{output_folder}/{prefix}_tDNAs.csv")
    cds_df.to_csv(f"{output_folder}/{prefix}_core_CDSs.csv")
    panaroo_roary.to_csv(f"{output_folder}/{prefix}_coregenome_detail.csv")
  
