import subprocess
import shlex
import ast
import pandas as pd
from Bio import SeqIO


def read_groups_tsv(tsv_file):
    tuples = []
    with open(tsv_file, "r") as table:
        for line in table:
            tuples.append(tuple(line.strip().split("\t")))
    clusts_dict = {}
    for i in [j[0] for j in tuples]:
        clusts_dict[i] = []
    for k in tuples:
        clusts_dict[k[0]].append(k[1])
    return clusts_dict


def read_coord_file(coordinates_file, dict_clusters):
    all_lcbs = []
    # Read the coordinates from the file
    with open(coordinates_file, "r") as file:
        next(file)
        # Assuming each line in the file contains start and end coordinates separated by a tab or space
        for line in file:
            tpl = line.strip().split("\t")
            all_lcbs.append([(tpl[i], tpl[i + 1]) for i in range(0, len(tpl), 2)])
    # Initialize a list to store the conserved contexts
    conserved_lcbs = [i for i in all_lcbs if ("0", "0") not in i]
    conserved_lcbs = [i for i in conserved_lcbs if int(i[0][1]) - int(i[0][0]) > 0]
    dict_str_blocks = {}
    for index, value in enumerate(dict_clusters):
        dict_str_blocks[value] = [pair[index] for pair in conserved_lcbs]
    return dict_str_blocks


def check_overlap(conserved_blocks, dict_clusters, tdna_name):
    ur_dict = {}
    dr_dict = {}
    for value, ranges in conserved_blocks.items():
        target_start = int(value.split("_")[3])
        target_end = int(value.split("_")[4])
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


def check_coordinates(tdna_df):
    block_coords = list(zip([int(i[1]) for i in list(tdna_df.LCB_UP)], [int(i[2]) for i in list(tdna_df.LCB_DOWN)]))
    tdnas_coords = list(zip(list(tdna_df.start), list(tdna_df.end)))
    tdna_row_name = [row.name for index, row in tdna_df.iterrows()]
    block_corrected = {}
    for index, pair in enumerate(list(zip(block_coords, tdnas_coords))):
        if pair[0][0] - pair[1][0] < 500:
            ur_corrected = pair[1][0] - 500
        else:
            ur_corrected = pair[1][0]
        if pair[1][1] - pair[0][1] < 500:
            dr_corrected = pair[1][1] + 500
        else:
            dr_corrected = pair[1][1]
        block_corrected[tdna_row_name[index]] = (ur_corrected, dr_corrected)
    return block_corrected


def detect_df_variations(dict_clusters, tdna_name):
    mauve_output = f"{tdna_name}_lcbs"
    seqs_frags = [line.strip() for line in open(mauve_output, "r") if line[0] in [">", "="]]
    groups = [i.split(" ")[1].split(":") for i in "".join(seqs_frags).split("=") if
              len(i.split("> ")) <= 2 and len(i) > 0]
    dict_new_ranges = {}
    for group in groups:
        seq_id = list(dict_clusters.keys())[int(group[0]) - 1]
        tuple_coords = tuple([int(i) for i in group[1].split("-")])
        if tuple_coords[1] - tuple_coords[0] > 50:
            dict_new_ranges.update({seq_id: tuple_coords})
    return dict_new_ranges


def new_dr_lcbs(tdna_df, lcbs_df):
    print(tdna_df)
    DR_dfs_core = []
    DR_dfs_non_core = []
    DR_LCBs = list(tdna_df["DR_name"].unique())
    DR_LCBs_count = tdna_df["DR_name"].value_counts().to_dict()
    for lcb in DR_LCBs_count.keys():
        lcb_df = lcbs_df.loc[lcbs_df["ID"] == str(lcb)]
        if len(lcb_df) == len(tdna_df["strain"].unique()):
            DR_dfs_core.append(dict(zip(lcb_df.strain, lcb_df.DR_tuple)))
        else:
            for group_name, group_df in lcb_df.groupby("ID"):
                DR_dfs_non_core.append(dict(zip(group_df.strain, group_df.DR_tuple)))
    return DR_dfs_core, DR_dfs_non_core


def change_dr_down(row, dr_dict):
    if row.sense == "+":
        return dr_dict[row.strain]
    else:
        return row.LCB_DOWN


def change_dr_up(row, dr_dict):
    if row.sense == "-":
        return dr_dict[row.strain]
    else:
        return row.LCB_UP


def extract_blocks_coords(scheme_df, lcbs_df, out_dir, numcpu):
    UR_ranges = {}
    DR_ranges = {}
    for tdna_name, tdna_df in scheme_df.groupby("tdna_name"):
        sequences = []
        tdna_df = tdna_df.copy()
        if len(tdna_df) > 1:
            # Check DR names
            print(tdna_name)
            core_dr, non_core_dr = new_dr_lcbs(tdna_df, lcbs_df)
            print(core_dr, non_core_dr, "AHOY!")
            if len(core_dr[0]) == len(tdna_df["strain"].unique()):
                test_df["LCB_DOWN"] = test_df.apply(lambda row: change_dr_down(row, core_dr[0]), axis=1)
                test_df["LCB_UP"] = test_df.apply(lambda row: change_dr_up(row, core_dr[0]), axis=1)
            else:
                print("AHOY!")

            # Extract FASTAs
            block_coords = check_coordinates(tdna_df)
            for index, row in tdna_df.iterrows():
                # fasta kintun formated
                fasta_file_path = f"{out_dir}/{row.strain}_annot/{row.strain}.fasta"
                record = SeqIO.read(fasta_file_path, "fasta")
                block_start, block_end = block_coords[index]
                if row.sense == "+":
                    sequence = record.seq[block_start:block_end]
                    tdna_start = int(row.start) - block_start
                    tdna_end = int(row.end) - block_start
                else:
                    sequence = record.seq[block_start:block_end].reverse_complement()
                    tdna_start = block_end - int(row.end)
                    tdna_end = block_end - int(row.start)
                sequences.append(f">{tdna_name}_{index}_{tdna_start}_{tdna_end}\n{sequence}")
            fasta_out_path = f"{tdna_name}_slices.fasta"
            with open(fasta_out_path, "w") as fasta_out:
                fasta_out.write("\n".join(sequences) + "\n")

            # Run mmseqs
            cmd1 = f"mmseqs easy-cluster {tdna_name}_slices.fasta {tdna_name}_clust tmp -c 0.99 --cov-mode 0 --threads {numcpu}"
            sp1 = subprocess.run(cmd1, shell=True, stdout=subprocess.DEVNULL)
            dict_clusters = read_groups_tsv(f"{tdna_name}_clust_cluster.tsv")
            if len(dict_clusters.keys()) > 1:
                cmd2 = f"progressiveMauve --output={tdna_name}_lcbs {tdna_name}_clust_rep_seq.fasta"
                sp2 = subprocess.run(cmd2, shell=True, stdout=subprocess.DEVNULL)
                ur_blocks, dr_blocks = check_overlap(read_coord_file(f"{tdna_name}_lcbs.backbone", dict_clusters),
                                                     dict_clusters, tdna_name)
                UR_ranges.update({int(i.split("_")[2]): (
                int(i.split("_")[3]) - ur_blocks[value][1], ur_blocks[value][1] - ur_blocks[value][0]) for index, value
                                  in enumerate(dict_clusters.keys()) for i in dict_clusters[value]})
                DR_ranges.update({int(i.split("_")[2]): (
                dr_blocks[value][0] - int(i.split("_")[4]), dr_blocks[value][1] - dr_blocks[value][0]) for index, value
                                  in enumerate(dict_clusters.keys()) for i in dict_clusters[value]})
            else:
                UR_ranges.update(
                    {int(i.split("_")[2]): (1, 1000) for index, value in enumerate(dict_clusters.keys()) for i in
                     dict_clusters[value]})
                DR_ranges.update(
                    {int(i.split("_")[2]): (1, 1000) for index, value in enumerate(dict_clusters.keys()) for i in
                     dict_clusters[value]})
        else:
            for index, row in tdna_df.iterrows():
                UR_ranges.update({index: (1, 1000)})
                DR_ranges.update({index: (1, 1000)})
    return UR_ranges, DR_ranges


def solve_ur_cooordinates(row):
    tdna_start = row.start
    tdna_end = row.end
    ur_dist, ur_size = row.UR_dist_size
    if row.sense == "+":
        if ur_size < 500 or ur_size > 5000:
            ur_coordinates = (tdna_start - ur_dist - 500, tdna_start - ur_dist)
        else:
            ur_coordinates = (tdna_start - ur_dist - ur_size, tdna_start - ur_dist)
    else:
        if ur_size < 500 or ur_size > 5000:
            ur_coordinates = (tdna_end + ur_dist, tdna_end + ur_dist + 500)
        else:
            ur_coordinates = (tdna_end + ur_dist, tdna_end + ur_dist + 500)
    return ur_coordinates


def solve_dr_cooordinates(row):
    tdna_start = row.start
    tdna_end = row.end
    dr_dist, dr_size = row.DR_dist_size
    if row.sense == "+":
        if dr_size < 500 or dr_size > 5000:
            dr_coordinates = (tdna_end + dr_dist, tdna_end + dr_dist + 500)
        else:
            dr_coordinates = (tdna_end + dr_dist, tdna_end + dr_dist + dr_size)
    else:
        if dr_size < 500 or dr_size > 5000:
            dr_coordinates = (tdna_start - dr_dist - 500, tdna_start - dr_dist)
        else:
            dr_coordinates = (tdna_start - dr_dist - dr_size, tdna_start - dr_dist)
    return dr_coordinates


def apply_nom(row, dict_ctx):
    range_ctx = dict_ctx[row.name]
    return range_ctx

def predict_gis(out_dir, scheme_file, lcbs_file, prefix, numcpu):
    # Import scheme as a DataFrame
    lcbs_df = pd.read_csv(lcbs_file)
    lcbs_df["IDs"] = lcbs_df.apply(lambda x: x["ID"].split("_")[0], axis=1)
    lcbs_df["DR_tuple"] = lcbs_df.apply(lambda row: (row.ID, row.start, row.end), axis=1)

    scheme_df = pd.read_csv(scheme_file, index_col=0)
    scheme_df['LCB_UP'] = scheme_df['LCB_UP'].apply(ast.literal_eval)
    scheme_df['LCB_DOWN'] = scheme_df['LCB_DOWN'].apply(ast.literal_eval)


    # Extract blocks coordinates
    UR_coords, DR_coords = extract_blocks_coords(scheme_df, lcbs_df, out_dir, numcpu)
    import json
    with open("ur_dict", "w") as file:
        json.dump(UR_coords, file)
    with open("dr_dict", "w") as file:
        json.dump(DR_coords, file)


    # Apply nom to DataFrame columns
    scheme_df["UR_dist_size"] = scheme_df.apply(lambda row: apply_nom(row, UR_coords), axis=1)
    scheme_df["DR_dist_size"] = scheme_df.apply(lambda row: apply_nom(row, DR_coords), axis=1)

    # Solve coordinates
    scheme_df["UR_coordinates"] = scheme_df.apply(lambda row: solve_ur_cooordinates(row), axis=1)
    scheme_df["DR_coordinates"] = scheme_df.apply(lambda row: solve_dr_cooordinates(row), axis=1)