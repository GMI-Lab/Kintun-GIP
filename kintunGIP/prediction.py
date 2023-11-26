import subprocess
import shlex
import ast
import pandas as pd
from Bio import SeqIO
import glob2
import os
from tqdm import tqdm

def read_groups_tsv(tsv_file):
    """
    This funtion reads groups from the TSV file output of mmseqs2,
    then transform it in a dictionary
    """
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
    """
    This function reads coordinates files from the progressiveMauve output
    and organize it based on their coordinates
    """
    all_lcbs = []

    with open(coordinates_file, "r") as file:
        next(file)
        # Assuming each line in the file contains start and end coordinates separated by a tab or space
        for line in file:
            tpl = line.strip().replace("-", "").split("\t")
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


def check_coordinates(row, window, genome_size):
    block_coords = (int(row.LCB_UP[1]), int(row.LCB_DOWN[2]))
    tdnas_coords = (row.start, row.end)
    pair = (block_coords, tdnas_coords)
    if pair[1][0] - pair[0][0] < window:
        ur_corrected = pair[1][0] - window
    else:
        ur_corrected = pair[0][0]
    if pair[1][1] - pair[0][1] < window:
        dr_corrected = pair[1][1] + window
    else:
        dr_corrected = pair[0][1]
    return ur_corrected, dr_corrected


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
    DR_dfs_core = []
    DR_dfs_non_core = []
    DR_LCBs = list(tdna_df["DR_name"].unique())
    DR_LCBs_count = tdna_df["DR_name"].value_counts().to_dict()
    for lcb in DR_LCBs_count.keys():
        lcb_df = lcbs_df.loc[lcbs_df["IDs"] == str(lcb)]
        if len(lcb_df) >= len(tdna_df["strain"].unique()):
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
        tdna_df = tdna_df.copy()
        if len(tdna_df) > 1:
            # Check DR names
            core_dr, non_core_dr = new_dr_lcbs(tdna_df, lcbs_df)
            if len(core_dr) != 0 and len(core_dr[0]) >= len(tdna_df["strain"].unique()):
                tdna_df["LCB_DOWN"] = tdna_df.apply(lambda row: change_dr_down(row, core_dr[0]), axis=1)
                tdna_df["LCB_UP"] = tdna_df.apply(lambda row: change_dr_up(row, core_dr[0]), axis=1)
            else:
                print(f"ERROR with {tdna_name}, no conserved UR or DR detected")
            # First iteration
            window = 20000
            ur_blocks, dr_blocks = {}, {}
            pbar = tqdm(desc="tmDNAs surroundings refinement", total = 500000)
            while window < 500000:
                sequences = []
                for index, row in tdna_df.iterrows():
                    # fasta kintun formated
                    fasta_file_path = f"{out_dir}/{row.strain}_annot/{row.strain}.fasta"
                    record = SeqIO.read(fasta_file_path, "fasta")
                    block_start, block_end = check_coordinates(row, window, len(record.seq))
                    if (0 - block_start + len(record.seq)) % len(record.seq) < (
                            block_end - block_start + len(record.seq)) % len(record.seq) + 1:
                        if row.sense == "+":
                            sequence = record.seq[block_start:] + record.seq[:block_end]
                            tdna_start = (int(row.start) - block_start + len(record.seq)) % len(record.seq) + 1
                            tdna_end = (int(row.end) - block_start + len(record.seq)) % len(record.seq) + 1
                        else:
                            sequence = (record.seq[block_start:] + record.seq[:block_end]).reverse_complement()
                            tdna_start = (block_end - int(row.end) + len(record.seq)) % len(record.seq) + 1
                            tdna_end = (block_end - int(row.start) + len(record.seq)) % len(record.seq) + 1
                    else:
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
                cmd1 = f"mmseqs easy-linclust {tdna_name}_slices.fasta {tdna_name}_clust tmp -c 0.90 --cov-mode 0 --threads {numcpu}"
                sp1 = subprocess.run(cmd1, shell=True, stdout=subprocess.DEVNULL)
                dict_clusters = read_groups_tsv(f"{tdna_name}_clust_cluster.tsv")
                if len(dict_clusters.keys()) > 1:
                    cmd2 = f"progressiveMauve --collinear --output={tdna_name}_lcbs {tdna_name}_clust_rep_seq.fasta"
                    sp2 = subprocess.run(cmd2, shell=True, stdout=subprocess.DEVNULL)
                    ur_blocks, dr_blocks = check_overlap(read_coord_file(f"{tdna_name}_lcbs.backbone", dict_clusters),
                                                         dict_clusters, tdna_name)
                if len(ur_blocks) != len(dict_clusters) or len(dr_blocks) != len(dict_clusters):
                    window += 20000
                    pbar.update(20000)
                    continue
                else:
                    break
            # Check coordinates
            dict_clusters = read_groups_tsv(f"{tdna_name}_clust_cluster.tsv")
            if len(dict_clusters.keys()) > 1:
                UR_ranges.update({int(i.split("_")[2]): (
                    int(i.split("_")[3]) - ur_blocks[value][1], ur_blocks[value][1] - ur_blocks[value][0]) for
                    index, value
                    in enumerate(dict_clusters.keys()) for i in dict_clusters[value]})
                DR_ranges.update({int(i.split("_")[2]): (
                    dr_blocks[value][0] - int(i.split("_")[4]), dr_blocks[value][1] - dr_blocks[value][0]) for
                    index, value
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
        for f in glob2.glob(f"{tdna_name}_*"):
            os.remove(f)
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


def check_ending(start_position, end_position, chromosome_size):
    # Ensure the end position is within the circular range
    if start_position < 0:
        start_position = chromosome_size + start_position
    elif start_position > chromosome_size:
        start_position = start_position - chromosome_size
    else:
        start_position = start_position
    # Ensure the end position is within the circular range
    if end_position < 0:
        end_position = chromosome_size + end_position
    elif end_position > chromosome_size:
        end_position = end_position - chromosome_size
    else:
        end_position = end_position
    return start_position, end_position


def final_coordinates_check(scheme_df, out_dir):
    coords_dict = {}
    for strain, strain_df in scheme_df.groupby("strain"):
        fasta_file_path = f"{out_dir}/{strain}_annot/{strain}.fasta"
        record = SeqIO.read(fasta_file_path, "fasta")
        for index, row in strain_df.iterrows():
            # tmdna
            start, end = check_ending(int(row.start), int(row.end), len(record.seq))
            # UR
            ur_coordinates = row.UR_coordinates
            ur_start, ur_end = check_ending(ur_coordinates[0], ur_coordinates[1], len(record.seq))
            # DR
            dr_coordinates = row.DR_coordinates
            dr_start, dr_end = check_ending(dr_coordinates[0], dr_coordinates[1], len(record.seq))
            coords_dict[f"{row.tdna_name}-{row.strain}"] = (start, end, ur_start, ur_end, dr_start, dr_end)
    return coords_dict


def add_coordinates_columns(checked_coordinates, row):
    values = checked_coordinates[f"{row.tdna_name}-{row.strain}"]
    return pd.Series([values[2], values[3], values[4], values[5]],
                     index=["ur_start", "ur_end", "dr_start", "dr_end"])


def check_distance(row):
    ur_distance = abs(row.ur_end - row.start)
    dr_distance = abs(row.dr_start - row.end)

    ur_info = "Empty"
    dr_info = "Empty"

    if ur_distance > 3000:
        ur_info = (min([row.start, row.ur_end]), max([row.start, row.ur_end]), abs(row.ur_end - row.start))

    if dr_distance > 3000:
        dr_info = (min([row.end, row.dr_start]), max([row.end, row.dr_start]), abs(row.dr_start - row.end))

    return ur_info, dr_info


def predict_gis(out_dir, scheme_file, lcbs_file, numcpu):
    # Import scheme as a DataFrame
    lcbs_df = pd.read_csv(lcbs_file)
    lcbs_df["IDs"] = lcbs_df.apply(lambda x: x["ID"].split("Block_")[1], axis=1)
    lcbs_df["DR_tuple"] = lcbs_df.apply(lambda row: (row.ID, row.start, row.end), axis=1)

    scheme_df = pd.read_csv(scheme_file, index_col=0)
    scheme_df['LCB_UP'] = scheme_df['LCB_UP'].apply(ast.literal_eval)
    scheme_df['LCB_DOWN'] = scheme_df['LCB_DOWN'].apply(ast.literal_eval)

    # Extract blocks coordinates
    UR_coords, DR_coords = extract_blocks_coords(scheme_df, lcbs_df, out_dir, numcpu)

    # Apply nom to DataFrame columns
    scheme_df["UR_dist_size"] = scheme_df.apply(lambda row: apply_nom(row, UR_coords), axis=1)
    scheme_df["DR_dist_size"] = scheme_df.apply(lambda row: apply_nom(row, DR_coords), axis=1)

    # Solve coordinates
    scheme_df["UR_coordinates"] = scheme_df.apply(lambda row: solve_ur_cooordinates(row), axis=1)
    scheme_df["DR_coordinates"] = scheme_df.apply(lambda row: solve_dr_cooordinates(row), axis=1)

    checked_coordinates = final_coordinates_check(scheme_df, out_dir)
    scheme_df[["ur_start", "ur_end", "dr_start", "dr_end"]] = scheme_df.apply(lambda row: add_coordinates_columns(checked_coordinates, row), axis=1)

    columns_to_keep = ['strain', 'start', 'end', 'sense', 'tdna_name', 'ur_start', 'ur_end', 'dr_start', 'dr_end']
    scheme_df = scheme_df[columns_to_keep]

    scheme_df[['ur_condition', 'dr_condition']] = scheme_df.apply(check_distance, axis=1, result_type='expand')

    scheme_df.to_csv(f"{scheme_file.replace('_tdna_scheme.csv','')}_tmDNAS_PILRs.csv")
