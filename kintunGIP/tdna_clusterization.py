import glob
from Bio import SeqIO
import os
import pandas as pd
import networkx as nx
from networkx.algorithms.community import louvain_communities
import string
from collections import Counter
from sklearn.metrics import jaccard_score

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

def find_neighbors(tRNA_row, cds_df, gene_ids):
    # Filter CDS rows based on the same "Strand"
    neighbors_cds = cds_df[cds_df['File'] == tRNA_row['File']]

    # Filter CDS rows based on their position relative to tRNA
    up_neighbors = neighbors_cds[
        (neighbors_cds['End'] < tRNA_row['Start']) |
        ((neighbors_cds['Start'] > tRNA_row['End']) & (neighbors_cds['Start'] <= tRNA_row['Start']))
    ].tail(3)
    if len(up_neighbors) < 2:
        up_neighbors = neighbors_cds.tail(3)

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


def cds_parse_genbank_files_to_dataframe(genbank_files, core_tags):
    data = {
        'File': [],
        'ID': [],
        'Type': [],
        'Start': [],
        'End': [],
        'Strand': []
    }

    for genbank_file in genbank_files:
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


def jaccard_similarity(list1, list2):
    s1 = set(list1)
    s2 = set(list2)
    if (len(s1.union(s2))) != 0:
        return float(len(s1.intersection(s2)) / (len(s1.union(s2))))
    else:
        return 0.0

def calculate_weight_hamming(data1, data2):
    # Your custom logic to calculate Hamming similarity for each attribute
    # hamming_similarity = sum(a == b for a, b in zip(data1, data2)) / len(data1)
    similarity = jaccard_similarity(data1[0], data2[0])
    if data1[2] == data2[2]:
        similarity += 1.0
    return similarity

def create_dict_ctxs(df):
    dict_ctx = {}
    for group_name, df_group in df.groupby("ANT"):
        if len(df_group) > 1:
            # Create an undirected graph
            G = nx.Graph()

            # Add nodes with row ID as node and attributes UR_name, tdnas_neigh, and strain
            for index, row in df_group.iterrows():
                node_id = row.name
                if row.Strand == 1:
                    ur_name = row['UP_neigh']  # Convert to set for easy intersection
                    dr_name = row['DOWN_neigh']  # Convert to set for easy intersection
                else:
                    ur_name = row['DOWN_neigh']  # Convert to set for easy intersection
                    dr_name = row['UP_neigh']  # Convert to set for easy intersection
                tdnas_neigh = row['tdnas_neigh']
                strain = row['File']

                G.add_node(node_id, UR_name=ur_name, DR_name=dr_name, tdnas_neigh=tdnas_neigh, strain=strain)

            # Add edges between nodes with a weight based on conditions
            # Add edges between nodes with a weight based on conditions
            threshold = 1.2# Adjust the threshold as needed

            for node1, data1 in G.nodes(data=True):
                for node2, data2 in G.nodes(data=True):
                    if node1 != node2:
                        # Calculate weight based on conditions
                        weight = calculate_weight_hamming([data1["UR_name"], data1["DR_name"], data1["tdnas_neigh"]],
                                                          [data2["UR_name"], data2["DR_name"], data2["tdnas_neigh"]])

                        # Add edge only if the weight is above the threshold
                        if weight > threshold:
                            G.add_edge(node1, node2, weight=weight)

            res = 0.8
            while True:
                communities = list(louvain_communities(G, weight='weight', resolution=res))

                # Check if all communities have unique strain values
                unique_strains = set()
                has_duplicate_strains = False

                for comm in communities:
                    community_strains = set(
                        node_data['strain'] for node, node_data in G.subgraph(comm).nodes(data=True))
                    if any(strain in unique_strains for strain in community_strains):
                        has_duplicate_strains = True
                        break
                    unique_strains.update(community_strains)

                res = res + 0.1
                if not has_duplicate_strains or res > 1.2:
                    break

            list_comm = []

            for i, community_set in enumerate(communities, 1):
                # Check for repeating strains within the community
                strain_counts = {}
                for node, node_data in G.subgraph(community_set).nodes(data=True):
                    strain = node_data['strain']
                    strain_counts[strain] = strain_counts.get(strain, 0) + 1

                # Identify nodes with repeating strains
                nodes_to_reassign = [node for node, node_data in G.subgraph(community_set).nodes(data=True) if
                                     strain_counts[node_data['strain']] > 1]

                # If there are nodes with repeating strains, reassign them to new communities
                if nodes_to_reassign:
                    new_communities = list(
                        louvain_communities(G.subgraph(community_set).copy(), weight='weight', resolution=0.9))
                    list_comm = list_comm + new_communities
                    for j, new_community_set in enumerate(new_communities, 1):
                        # Update the community assignment in the original graph
                        for node in new_community_set:
                            G.nodes[node]['community'] = str(i) + "_" + str(
                                j)  # Using a unique identifier for reassigned communities

                else:
                    # If no repeating strains, continue with the original community assignment
                    list_comm = list_comm + [community_set]
                    for node in community_set:
                        G.nodes[node]['community'] = i

            list_comm_sorted = sorted(list_comm, key=len, reverse=True)

            community_mapping = {}
            for i, community in enumerate(list_comm_sorted):
                for node in community:
                    community_mapping[node] = i

            labels = [community_mapping[j] for j, row in df_group.iterrows()]

        else:
            labels = [0]

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

def tdna_clusterization(input_folder, output_folder, prefix, format):
    # Example: Parse a list of GenBank files and create a DataFrame
    genbank_files_list = glob.glob(f"{input_folder}/*.{format}")
    tdnas_df = parse_genbank_files_to_dataframe(genbank_files_list)

    panaroo_roary = pd.read_csv(f"{input_folder}/gene_presence_absence_roary.csv", header=0, low_memory=False)
    panaroo_roary = panaroo_roary.loc[(panaroo_roary["No. isolates"] > 0.9 * len(genbank_files_list)) & (
                panaroo_roary["No. sequences"] > 0.9 * len(genbank_files_list))]

    # Create a dictionary
    core_tags = {}

    # Iterate over columns starting from the 15th column
    for column in panaroo_roary.columns[14:]:
        values = panaroo_roary[column].apply(lambda x: x.split(';') if ';' in str(x) else [x]).tolist()
        core_tags[column] = values

    gene_ids = {}
    for index, row in panaroo_roary.iterrows():
        key = row['Gene']
        values = {}

        # Iterate over columns starting from the 15th column
        for column, value in row.iloc[14:].items():
            # Split values if semicolon is present
            value_list = str(value).split(';')
            for split_value in value_list:
                gene_ids[split_value] = key

    cds_df = cds_parse_genbank_files_to_dataframe(genbank_files_list, core_tags)
    tdnas_df["ANT"] = tdnas_df["Product"].apply(apply_anticodon_mapping)
    tdnas_df = tdnas_df[tdnas_df['Product'] != 'tRNA-Xxx']
    neighbors_series = tdnas_df.apply(lambda row: find_neighbors(row, cds_df, gene_ids), axis=1)
    tdnas_df['UP_neigh'] = neighbors_series.apply(lambda x: x[0])
    tdnas_df['DOWN_neigh'] = neighbors_series.apply(lambda x: x[1])
    tdnas_df['UP_neigh_ori'] = neighbors_series.apply(lambda x: x[2])
    tdnas_df['DOWN_neigh_ori'] = neighbors_series.apply(lambda x: x[3])
    tdnas_df["tdnas_neigh"] = tdnas_df.apply(tdnas_neigh_finder, args=(tdnas_df,), axis=1)
    tdnas_df.to_csv(f"{output_folder}/{prefix}_prev_tDNAs.csv")
    dict_ctx = create_dict_ctxs(tdnas_df)
    tdnas_df["tdna_name"] = tdnas_df.apply(apply_nom, args=(dict_ctx, prefix,), axis=1)
    tdnas_df.to_csv(f"{output_folder}/{prefix}_tDNAs.csv")
    cds_df.to_csv(f"{output_folder}/{prefix}_core_CDSs.csv")
    panaroo_roary.to_csv(f"{output_folder}/{prefix}_coregenome_detail.csv")
