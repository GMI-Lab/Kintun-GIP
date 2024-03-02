import pandas as pd
from ast import literal_eval
import os
from Bio import SeqIO
import subprocess


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


def check_cds(tdnas_df, cds_df, core_df):
    UR_left = {}
    DR_right = {}
    for group_name, group_df in tdnas_df.groupby("tdna_name"):
        UR_groups = []
        DR_groups = []
        for index, row in group_df.iterrows():
            if row.Strand == 1:
                for i in literal_eval(row.UP_neigh):
                    UR_groups.append(i.replace("_Yes", "").replace("_No", ""))
                for j in literal_eval(row.DOWN_neigh):
                    DR_groups.append(j.replace("_Yes", "").replace("_No", ""))
            else:
                for i in literal_eval(row.DOWN_neigh):
                    UR_groups.append(i.replace("_Yes", "").replace("_No", ""))
                for j in literal_eval(row.UP_neigh):
                    DR_groups.append(j.replace("_Yes", "").replace("_No", ""))
        UR_groups = list(set(UR_groups))
        DR_groups = list(set(DR_groups))

        for index, row in group_df.iterrows():
            UR_coords = []
            DR_coords = []
            for k in UR_groups:
                condition_row = cds_df['ID'] == core_df.at[k, row.File]
                result = [tuple(x) for x in cds_df.loc[condition_row, ['Start', 'End']].to_numpy()]
                if len(result) == 1:
                    UR_coords.append(result[0])
            for l in DR_groups:
                condition_row = cds_df['ID'] == core_df.at[l, row.File]
                result = [tuple(x) for x in cds_df.loc[condition_row, ['Start', 'End']].to_numpy()]
                if len(result) == 1:
                    DR_coords.append(result[0])
            if row.Strand == 1:
                UR_left[row.name] = min([range1 for range1 in UR_coords if
                                         check_if_overlap_circular(range1, (row.Start - 250000, row.Start),
                                                                   row.Chr_size)])
                DR_right[row.name] = max([range1 for range1 in DR_coords if
                                          check_if_overlap_circular(range1, (row.End, row.End + 250000), row.Chr_size)])
            else:
                UR_left[row.name] = min([range1 for range1 in DR_coords if
                                         check_if_overlap_circular(range1, (row.Start - 250000, row.Start),
                                                                   row.Chr_size)])
                DR_right[row.name] = max([range1 for range1 in UR_coords if
                                          check_if_overlap_circular(range1, (row.End, row.End + 250000), row.Chr_size)])
    return UR_left, DR_right


def apply_nom(row, dict_ctx):
    tuple_value = dict_ctx[row.name]
    return tuple_value


def extract_sequence(genbank_file, strand, start, end, chromosome_size):
    records = SeqIO.parse(genbank_file, "genbank")

    for record in records:
        # Extract the circular sequence
        circular_sequence = record.seq + record.seq

        # Extract the region
        if end < start:
            end = end + chromosome_size

        if strand == 1:
            region_sequence = circular_sequence[start:end]
        else:
            region_sequence = circular_sequence[start:end].reverse_complement()
        return region_sequence


def create_dict_row(block, dct_codes):
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
            "strain": dct_codes[lcb[0]],
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


def create_lcbs_df(blocks_file):
    df = pd.DataFrame(columns=["strain", "method", "type", "start", "end", "na1", "sense", "na2", "ID"])
    with open(blocks_file, "r") as lcbs_file:
        fields = "".join(lcbs_file.readlines()).split(f"{80 * '-'}\n")
        seqs = [i.strip().split("\n") for i in fields if not "Block #" in i][0][1:]
        groups = [i.strip().split("\n") for i in fields if "Block #" in i]
        conserved_lcb = [i for i in groups]
    dct_codes = {}
    for i in seqs:
        line = i.strip().split("\t")
        dct_codes[line[0]] = line[-1]
    rows_lcb = [create_dict_row(block, dct_codes) for block in conserved_lcb]
    new_rows = [row for group in rows_lcb for row in group]
    df_new = pd.DataFrame(new_rows, columns=["strain", "method", "type", "start", "end", "na1", "sense", "na2", "ID"])
    df = pd.concat([df, df_new])
    for group_name, df_group in df.groupby("strain"):
        df_group['ID'] = df_group['ID'] + '_' + df_group.groupby('ID').cumcount().add(1).astype(str)
        df.update(df_group)
    return df


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


def delete_non_core(lcbdf):
    groups = lcbdf.groupby("ID", group_keys=False)
    filtered_groups = groups.filter(lambda x: len(x) == lcbdf['strain'].nunique())
    new_df = filtered_groups.groupby("ID", group_keys=False).apply(lambda x: pd.concat([x]))
    new_df.rename(columns={'ID': 'IDs'}, inplace=True)
    return new_df


def adjust_ranges(row):
    target_start, target_end = (int(row.strain.split("-")[2]), int(row.strain.split("-")[3]))
    lcb_start, lcb_end = (row.start, row.end)
    # Check for overlap
    if lcb_end in range(target_start, target_end) and lcb_start <= target_start:
        lcb_start = lcb_start
        lcb_end = target_start - 1
        return (lcb_start, lcb_end)
    elif lcb_start in range(target_start, target_end) and lcb_end >= target_end:
        lcb_start = target_end + 1
        lcb_end = lcb_end
        return (lcb_start, lcb_end)
    elif lcb_start <= target_start and lcb_end >= target_end:
        part1 = (lcb_start, target_start - 1)
        part2 = (target_end + 1, lcb_end)
        return (part1, part2)
    else:
        return (lcb_start, lcb_end)


def is_nested_tuple(input_tuple):
    # Check if each element in the tuple is an instance of tuple
    return isinstance(input_tuple, tuple) and all(isinstance(item, tuple) for item in input_tuple)


def find_nearest_lower_range(target_range, ranges_list):
    lower_ranges = [rng for rng in ranges_list if rng[0] < target_range[0]]
    if not lower_ranges:
        return "None"  # No lower ranges found

    nearest_lower_range = max(lower_ranges, key=lambda rng: rng[0])
    return nearest_lower_range


def find_nearest_upper_range(target_range, ranges_list):
    upper_ranges = [rng for rng in ranges_list if rng[0] > target_range[0]]
    if not upper_ranges:
        return "None"  # No upper ranges found

    nearest_upper_range = min(upper_ranges, key=lambda rng: rng[0])
    return nearest_upper_range


def conserved_ctx_def(row, lcbs_df, group_name):
    if row.Strand == 1:
        target_start = row.Start - row.bor_left[0]
        target_end = row.End - row.bor_left[0]
    else:
        target_start = row.bor_right[1] - row.End
        target_end = row.bor_right[1] - row.Start
    lcbs_cons = lcbs_df.loc[lcbs_df["strain"] == f"{group_name}-{row.name}-{target_start + 1}-{target_end}"]
    CTX_tuples = [row.cor_coords for index, row in lcbs_cons.iterrows()]
    CTX_ranges = []
    for value in CTX_tuples:
        if is_nested_tuple(value):
            for i in value:
                if row.Strand == 1:
                    CTX_ranges.append((i[0] + row.bor_left[0], i[1] + row.bor_left[0]))
                else:
                    CTX_ranges.append((row.bor_right[1] - i[0], row.bor_right[1] - i[1]))
        else:
            if row.Strand == 1:
                CTX_ranges.append((value[0] + row.bor_left[0], value[1] + row.bor_left[0]))
            else:
                CTX_ranges.append((row.bor_right[1] - value[0], row.bor_right[1] - value[1]))
    UR_range = find_nearest_lower_range((row.Start, row.End), CTX_ranges)
    DR_range = find_nearest_upper_range((row.Start, row.End), CTX_ranges)
    return UR_range, DR_range


def check_gaps(target_range, other_ranges, threshold=3000):
    target_start, target_end = sorted(target_range)
    other_ranges = list(other_ranges)
    if len(other_ranges[0]) == 2 or len(other_ranges[1]) == 2:
        if len(other_ranges[0]) == 2 and len(other_ranges[1]) == 2:
            # Check for gaps on the left
            left_start, left_end = other_ranges[0]
            if target_start - left_end > threshold:
                left_gap = (left_end + 1, target_start - 1)
            else:
                left_gap = "None"
            # Check the gap between the target end and the last range
            right_start, right_end = other_ranges[1]
            if right_start - target_end > threshold:
                right_gap = (target_end + 1, right_start - 1)
            else:
                right_gap = "None"
            return left_gap, right_gap
        # 3
        elif len(other_ranges[0]) != 2 and len(other_ranges[1]) == 2:
            right_start, right_end = other_ranges[1]
            if right_start - target_end > threshold:
                right_gap = (target_end + 1, right_end)
            else:
                right_gap = "None"
            return "None", right_gap
        # 4
        elif len(other_ranges[0]) == 2 and len(other_ranges[1]) != 2:
            left_start, left_end = other_ranges[0]
            if target_start - left_start > 2000:
                left_gap = (left_end + 1, target_start)
            else:
                left_gap = "None"
            return left_gap, "None"
        # Return the results
    else:
        return "None", "None"


def pre_predict_PILR(tdnas_df, out_dir, input_dir):
    ur_dr_dict = {}
    for group_name, group_df in tdnas_df.groupby("tdna_name"):
        if len(group_df) >= 2:
            fasta_ctxs = f"{out_dir}/{group_name}_ctx.fa"
            if not os.path.isfile(fasta_ctxs):
                with open(fasta_ctxs, "w") as out_file:
                    for index, row in group_df.iterrows():
                        if row.Strand == 1:
                            target_coords = (abs(row.Start - row.bor_left[0]), abs(row.End - row.bor_left[0]))
                            target_start = min(target_coords)
                            target_end = max(target_coords)
                        else:
                            target_coords = (abs(row.bor_right[1] - row.End), abs(row.bor_right[1] - row.Start))
                            target_start = min(target_coords)
                            target_end = max(target_coords)
                        out_file.write(
                            f">{group_name}-{row.name}-{target_start + 1}-{target_end}\n{str(extract_sequence(f'{input_dir}/{row.File}.gbff', row.Strand, row.bor_left[0], row.bor_right[1], row.Chr_size))}\n")
            if not os.path.isfile(f"{out_dir}/sibelia_{group_name}/300/blocks_coords.txt"):
                print(group_name, "SibeliaZ")
                cmd2 = f"sibeliaz -k 11 -n -t 8 -o {out_dir}/sibelia_{group_name}/ {fasta_ctxs} ; maf2synteny -b 300 -o {out_dir}/sibelia_{group_name}/ {out_dir}/sibelia_{group_name}/blocks_coords.gff"
                sp2 = subprocess.run(cmd2, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            lcbs_df = create_lcbs_df(f"{out_dir}/sibelia_{group_name}/300/blocks_coords.txt")
            lcbs_df = delete_non_core(lcbs_df)
            lcbs_df["cor_coords"] = lcbs_df.apply(adjust_ranges, axis=1)
            for index, row in group_df.iterrows():
                ur_dr_dict[row.name] = conserved_ctx_def(row, lcbs_df, group_name)
        else:
            for index, row in group_df.iterrows():
                ur_dr_dict[row.name] = ((row.Start - 2001, row.Start - 1), (row.End + 1, row.End + 2001))
    return ur_dr_dict


def PILR_prediction(row, ur_dr_dict):
    target_range = (row.Start, row.End)
    locus_ranges = ur_dr_dict[row.name]
    left_gap, right_gap = check_gaps(target_range, locus_ranges)
    return (left_gap, right_gap)


def find_overlapping_sublists(ranges):
    overlapping_sublists = []
    current_sublist = [ranges[0]]
    for current_range in ranges[1:]:
        last_range = current_sublist[-1]
        if current_range[0] <= last_range[1]:  # Check for overlap
            current_sublist.append(current_range)
        else:
            if len(current_sublist) > 1:
                overlapping_sublists.append(current_sublist)
            current_sublist = [current_range]

    if len(current_sublist) > 1:
        overlapping_sublists.append(current_sublist)

    return overlapping_sublists


def get_intersection_range(overlapping_sublist):
    # Sort the ranges based on their start values
    sorted_ranges = sorted(overlapping_sublist, key=lambda x: x[0])
    # Initialize the intersection range
    intersection_start, intersection_end = sorted_ranges[0]
    # Iterate through the sorted ranges to find the intersection
    for start, end in sorted_ranges[1:]:
        if start > intersection_end:  # No overlap
            return None
        else:
            intersection_start = max(intersection_start, start)
            intersection_end = min(intersection_end, end)
    return (intersection_start, intersection_end)


def find_best_gene(intersection, gene_options):
    # Sort gene options based on the distance from the end coordinate of the intersection
    sorted_genes = sorted(gene_options, key=lambda gene: abs(gene[2] - intersection[1]))
    # The best gene is the one with the smallest distance to the end coordinate of the intersection
    best_gene = sorted_genes[0]
    return best_gene


def filt_overlap_PILRs(df):
    # First_left
    dict_coords_left = {}
    for group_name, group_df in df.groupby("File"):
        tdnas_names = [(row.name, int(row.Start), int(row.End), row.Strand) for index, row in group_df.iterrows() if
                       row.PILRs[0] != "None"]
        PILRs_left = [row.PILRs[0] for index, row in group_df.iterrows() if row.PILRs[0] != "None"]
        to_check = find_overlapping_sublists(PILRs_left)
        if len(to_check) > 0:
            for overlapping_sublist in to_check:
                intersection_range = get_intersection_range(overlapping_sublist)
                all_genes = [tdnas_names[PILRs_left.index(j)] for j in overlapping_sublist]
                best_gene = find_best_gene(intersection_range, all_genes)
                dict_coords_left[best_gene[0]] = intersection_range
                for gene in all_genes:
                    if gene[0] not in dict_coords_left.keys():
                        dict_coords_left[gene[0]] = "None"

                        # Second right
    dict_coords_right = {}
    for group_name, group_df in df.groupby("File"):
        tdnas_names = [(row.name, int(row.Start), int(row.End), row.Strand) for index, row in group_df.iterrows() if
                       row.PILRs[1] != "None"]
        PILRs_right = [row.PILRs[1] for index, row in group_df.iterrows() if row.PILRs[1] != "None"]
        to_check = find_overlapping_sublists(PILRs_right)
        if len(to_check) > 0:
            for overlapping_sublist in to_check:
                intersection_range = get_intersection_range(overlapping_sublist)
                all_genes = [tdnas_names[PILRs_right.index(j)] for j in overlapping_sublist]
                best_gene = find_best_gene(intersection_range, all_genes)
                dict_coords_right[best_gene[0]] = intersection_range
                for gene in all_genes:
                    if gene[0] not in dict_coords_right.keys():
                        dict_coords_right[gene[0]] = "None"

    return dict_coords_left, dict_coords_right


def check_overlap_PILR(row, dict_pilrs, num):
    if row.name in dict_pilrs.keys():
        return dict_pilrs[row.name]
    else:
        return row.PILRs[num]


def predict_gis(output_dir, input_dir, prefix):
    tdnas_df = pd.read_csv(f"{output_dir}/{prefix}_tDNAs.csv", header=0, index_col=0)
    cds_df = pd.read_csv(f"{output_dir}/{prefix}_core_CDSs.csv", header=0, index_col=0)
    core_df = pd.read_csv(f"{output_dir}/{prefix}_coregenome_detail.csv", header=0, index_col=1)
    left_dict, right_dict = check_cds(tdnas_df, cds_df, core_df)
    tdnas_df["bor_left"] = tdnas_df.apply(apply_nom, args=(left_dict,), axis=1)
    tdnas_df["bor_right"] = tdnas_df.apply(apply_nom, args=(right_dict,), axis=1)
    ur_dr_dict = pre_predict_PILR(tdnas_df, output_dir, input_dir)
    tdnas_df["PILRs"] = tdnas_df.apply(PILR_prediction, args=(ur_dr_dict,), axis=1)
    pilrs_cor_left, pilrs_cor_right = filt_overlap_PILRs(tdnas_df)
    tdnas_df["PILR_left"] = tdnas_df.apply(check_overlap_PILR, args=(pilrs_cor_left, 0,), axis=1)
    tdnas_df["PILR_right"] = tdnas_df.apply(check_overlap_PILR, args=(pilrs_cor_right, 0,), axis=1)
    tdnas_df.to_csv(f"{output_dir}/{prefix}_PILRs.csv")
