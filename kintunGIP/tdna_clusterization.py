# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 18:34:12 2020

@author: camilo
"""

# Imports
import pandas as pd
import os
import glob
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import FeatureLocation
from collections import Counter
import numpy as np
import logging
import ntpath
from tqdm.auto import tqdm
from sklearn.cluster import DBSCAN
import scipy.spatial
from collections import Counter
import string

# Dictionary with anticodon codes
dct_ant = dict(GLYACC='gly4', HISATG='his2', CYSACA='cys2', ARGACG='arg1', ASPATC='asp2',
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

# Functions
def apply_ant_dict(ID_info):
    tmdna_id = ""
    if "Name=tmRNA" in ID_info:
        tmdna_id = "ssrA"
    else:
        anticodon = ID_info.split("Name=")[1].replace("-", "")
        tmdna_id = dct_ant[anticodon.upper()]
    return tmdna_id


def accession_dct(list_files):
  str_dct = {}
  for in_file in list_files:
    with open(in_file,"r") as strain_file:
      for record in SeqIO.parse(in_file, "fasta"):
        str_dct[record.id] = ntpath.basename(in_file).replace(f".fasta", "")
  return str_dct

def check_core(no_strains, lcb_group):
  codes = [int(i.split("\t")[0]) for i in lcb_group[2:]]
  for num in range(1,no_strains+1):
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
            "type"  : "LCB-Core",
            "start" : lcb_start,
            "end"   : lcb_end,
            "na1"   : ".",
            "sense" : lcb[-4],
            "na2"   : ".",
            "ID"    : f"Name=Block_#{ident}"
        }
        new_rows.append(new_row)
    return new_rows


def create_lcbs_df(blocks_file, list_files):
    df = pd.DataFrame(columns=["strain","method","type","start","end","na1","sense","na2","ID"])
    dct_files = accession_dct(list_files)
    with open(blocks_file,"r") as lcbs_file:
        fields = "".join(lcbs_file.readlines()).split(f"{80*'-'}\n")
        # Extract information
        seqs = [i.strip().split("\n") for i in fields if not "Block #" in i][0][1:]
        groups = [ i.strip().split("\n") for i in fields if "Block #" in i]
        conserved_lcb = [i for i in groups if check_core(len(seqs),i)]
        # Obtain coords
    dct_codes = {}
    for i in seqs:
        line = i.strip().split("\t")
        dct_codes[line[0]] = line[-1]
    rows_lcb = [create_dict_row(block, dct_codes, dct_files) for block in conserved_lcb]
    new_rows = [row for group in rows_lcb for row in group]
    df = df.append(new_rows, ignore_index=True)

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
    groups = lcbdf.groupby("ID")
    filtered_groups = groups.filter(lambda x: len(x) >= lcbdf['strain'].nunique())
    new_df = filtered_groups.groupby("ID").apply(lambda x: pd.concat([x]))
    new_df.rename(columns = {'ID':'IDs'}, inplace = True)
    return new_df

  
def run_sibeliaz(list_files, output_folder, threads):
    cmd_sibeliaz = f"sibeliaz -k 11 -n -t {threads} -o {output_folder}/all_chr_sibelia/ {' '.join(list_files)} ; maf2synteny -b 50 -o {output_folder}/all_chr_sibelia/ {output_folder}/all_chr_sibelia/blocks_coords.gff"
    p = subprocess.run(cmd_sibeliaz, shell=True)

    
def seqlen_dct(list_files, file_ext):
  str_dct = {}
  for in_file in list_files:
    with open(in_file,"r") as strain_file:
      for record in SeqIO.parse(in_file, "fasta"):
        str_dct[ntpath.basename(in_file).replace(f".{file_ext}", "")] = len(record.seq)
  return str_dct


def overlap_at_start(r1, r2):
    """Check if two ranges overlap at the start."""
    return r1[0] == r2[0] or (r1[0] < r2[0] and r1[1] > r2[0])

  
def overlap_at_end(r1, r2):
    """Check if two ranges overlap at the end."""
    return r1[1] == r2[1] or (r1[0] < r2[1] and r1[1] > r2[1])

  
def check_if_overlap(range1, range2):
    return not (range1[1] < range2[0] or range2[1] < range1[0])

  
def find_nearest_range(point, ranges):
    """Find the nearest range in a list of ranges to a given point."""
    nearest_range = min(ranges, key=lambda r: abs(point - r[0]))
    return nearest_range

  
def up_lcb_finder(row, dct_len, lcbdf):
    len_seq = dct_len[row.strain]
    strain_lcb = lcbdf[lcbdf['strain'] == row.strain]
    strain_lcb = strain_lcb.sort_values(by=["start"])
    lcb_tuples_ori = [(start, end) for start,end in zip(strain_lcb['start'], strain_lcb['end'])]
    lcb_tuples =     [(start, end) for start,end in zip(strain_lcb['start'], strain_lcb['end'])]
    lcb_tuples.append((lcb_tuples[-1][0]-len_seq,lcb_tuples[-1][1]-len_seq))
    lcb_tuples.append((lcb_tuples[0][0]+len_seq,lcb_tuples[0][1]+len_seq))
    non_overlap = [lcb for lcb in lcb_tuples if not check_if_overlap(lcb, (row.start,row.end))]
    nearest_ur = find_nearest_range(row.start, [ranges for ranges in non_overlap if ranges[0] < row.start])
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
    lcb_tuples_ori = [(start, end) for start,end in zip(strain_lcb['start'], strain_lcb['end'])]
    lcb_tuples = [(start, end) for start,end in zip(strain_lcb['start'], strain_lcb['end'])]
    lcb_tuples.append((lcb_tuples[-1][0]-len_seq,lcb_tuples[-1][1]-len_seq))
    lcb_tuples.append((lcb_tuples[0][0]+len_seq,lcb_tuples[0][1]+len_seq))
    non_overlap = [lcb for lcb in lcb_tuples if not check_if_overlap(lcb, (row.start,row.end))]
    nearest_ur = find_nearest_range(row.start, [ranges for ranges in non_overlap if ranges[0] > row.end])
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
        return row["LCB_UP"][0].replace("Name=Block_","")
    else:
        return row["LCB_DOWN"][0].replace("Name=Block_","")

      
def DR_block_name(row):
    if row["sense"] == "+":
        return row["LCB_DOWN"][0].replace("Name=Block_","")
    else:
        return row["LCB_UP"][0].replace("Name=Block_","")

      
def check_if_overlap(tuple1,tuple2):
    x = range(tuple1[0],tuple1[1])
    y = range(tuple2[0],tuple2[1])
    if len(range(max(x[0], y[0]), min(x[-1], y[-1])+1)) > 0:
        return True
    else:
        return False

      
def tdnas_neigh_finder(row, dct_len, tmrnadf):
    len_seq = dct_len[row.strain]
    if row.sense == "+":
        range_neigh = (row.start-1000, row.start-1)
    else:
        range_neigh = (row.end+1, row.end+1000)
    strain_tmrnas = tmrnadf[tmrnadf['strain'] == row.strain]
    strain_tmrnas = strain_tmrnas.sort_values(by=["start"])
    tmrnas_tuples = [(start, end) for start,end in zip(strain_tmrnas['start'], strain_tmrnas['end'])]
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
        return LCB_UR[0]+"_YES"
    else:
        return LCB_UR[0]+"_NO"


def create_dict_ctxs(df):
    dict_ctx = {}
    for group_name, df_group in df.groupby("ANT"):
        if len(df_group) > 1:
            df_clust = df_group[["UR_name", "tdnas_neigh", "coin_sense"]]
            df_clust_binary = pd.get_dummies(df_clust)
            jaccard = scipy.spatial.distance.cdist(df_clust_binary, df_clust_binary, metric='jaccard')
            rows_distance = pd.DataFrame(jaccard, columns=df_clust_binary.index.values,
                                         index=df_clust_binary.index.values)
            # Instantiate the DBSCAN object with the desired hyperparameters
            dbscan = DBSCAN(eps=0.4, min_samples=2, metric='precomputed')
            # Perform DBSCAN clustering on the distance matrix
            labels = dbscan.fit_predict(rows_distance)
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

  
def modify_duplicate(df):
    # Group the DataFrame by "City" and apply a cumulative count
    counts = df.groupby(["strain", "tdna_name"]).cumcount() + 1
    # Create a boolean mask to only modify duplicated values
    mask = df.duplicated(subset=["strain", "tdna_name"], keep=False)
    # Concatenate the "Name" column with the counts, separated by an underscore, only for duplicated values
    df.loc[mask, "tdna_name"] = df.loc[mask, "tdna_name"] + "_" + counts.astype(str)
    return df

  
def create_feature(row):
    salida_features = []
    #Create tmDNA feature                       
    if row["start"] < row["end"]:
        start_pos = ExactPosition(row["start"])
        end_pos = ExactPosition(row["end"])
    else:
        start_pos = ExactPosition(row["end"])
        end_pos = ExactPosition(row["start"])
    feature_type = "tRNA"
    if row["sense"] == "+":
      sense = +1
    else:
      sense = -1
    feature_location = FeatureLocation(start_pos, end_pos, strand=sense)
    id_string = f'{row["tdna_name"]}'
    feature = SeqFeature(
    feature_location,
    type = feature_type,
    id = id_string,
    qualifiers={
        "locus_tag" : id_string,
        }
    )
    salida_features.append(feature)
    #Create UR feature
    #(Name=Block_#1871, 3280981, 3281075, -)
    ur_feat = row.LCB_UP
    if ur_feat[1] < ur_feat[2]:
        start_pos = ExactPosition(ur_feat[1])
        end_pos = ExactPosition(ur_feat[2])
    else:
        start_pos = ExactPosition(ur_feat[2])
        end_pos = ExactPosition(ur_feat[1])
    feature_type = "UR_BLOCK"
    if ur_feat[3] == "+":
      sense = +1
    else:
      sense = -1
    feature2_location = FeatureLocation(start_pos, end_pos, strand=sense)
    id_string = f'{ur_feat[0]}'
    feature2 = SeqFeature(
    feature2_location,
    type = feature_type,
    id = id_string,
    qualifiers={
        "locus_tag" : ur_feat[0] + "_" + row.tdna_name,
        }
    )
    salida_features.append(feature2)
    #Create DR feature
    dr_feat = row.LCB_DOWN
    if dr_feat[1] < dr_feat[2]:
        start_pos = ExactPosition(dr_feat[1])
        end_pos = ExactPosition(dr_feat[2])
    else:
        start_pos = ExactPosition(dr_feat[2])
        end_pos = ExactPosition(dr_feat[1])
    feature_type = "DR_BLOCK"
    if dr_feat[3] == "+":
      sense = +1
    else:
      sense = -1
    feature3_location = FeatureLocation(start_pos, end_pos, strand=sense)
    id_string = f'{dr_feat[0]}'
    feature3 = SeqFeature(
    feature3_location,
    type = feature_type,
    id = id_string,
    qualifiers={
        "locus_tag" : dr_feat[0] + "_" + row.tdna_name,
        }
    )
    salida_features.append(feature3)
    return salida_features

  
def create_genbanks(list_files, tmrnadf, file_ext):
    for fasta in tqdm(list_files):
        with open(fasta,"r") as fasta_in:
            for record in SeqIO.parse(fasta_in,"fasta"):
                new_features = []
                for index, row in tmrnadf.loc[tmrnadf['strain'] == ntpath.basename(fasta).replace(f".fasta", "")].iterrows():
                    new_features.append(create_feature(row))
                for i in new_features:
                    for j in i:
                        record.features.append(j)
                record.annotations["molecule_type"] = "DNA"
                record.annotations["form"] = "double stranded"
                record.annotations["topology"] = "circular"
                with open(fasta.replace(f".{file_ext}","_KintunClust.gb"),"w") as salida:
                    SeqIO.write(record,salida,"gb")

                    
def tdna_clusterization(input_folder, output_folder, file_ext, nom_ext, threads):
    # Create list of files .aragorn
    # list_files = glob.glob(f'{output_folder}/*/*.tmdnas', recursive=True)
    list_files_trnas = glob.glob(f'{output_folder}/*/*.tmdnas', recursive=True)
    list_files_fasta = glob.glob(f'{output_folder}/*/*.fasta', recursive=True)
    # Create tmdnas dataframe
    tmrnadf = create_tmdnas_df(list_files_trnas)
    # Run SibeliaZ with all genomes
    run_sibeliaz(list_files_fasta, output_folder, threads)
    # Create dataframe with LCBs
    lcbdf = create_lcbs_df(f"{output_folder}/all_chr_sibelia/50/blocks_coords.txt", list_files_fasta)
    lcbdf = delete_non_core(lcbdf)
    # Create chr length dict
    dct_len = seqlen_dct(list_files_fasta, file_ext)
    # Create DR and UR info columns
    tqdm.pandas(desc="Processing <--(t(m)DNA)")
    tmrnadf["LCB_UP"] = tmrnadf.progress_apply(lambda row: up_lcb_finder(row, dct_len, lcbdf), axis=1)
    tqdm.pandas(desc="Processing (t(m)DNA)-->")
    tmrnadf["LCB_DOWN"] = tmrnadf.progress_apply(lambda row: down_lcb_finder(row, dct_len, lcbdf), axis=1)
    # Correct UR and DR names by sense
    tqdm.pandas(desc="Processing...")
    tmrnadf["UR_name"] = tmrnadf.progress_apply(lambda row: UR_block_name(row), axis=1)
    tmrnadf["DR_name"] = tmrnadf.progress_apply(lambda row: DR_block_name(row), axis=1)
    tmrnadf["ANT"] = tmrnadf.progress_apply(lambda row: apply_ant_dict(row["ID"]), axis=1)
    tmrnadf["tdnas_neigh"] = tmrnadf.progress_apply(lambda row: tdnas_neigh_finder(row, dct_len, tmrnadf), axis=1)
    tmrnadf["coin_sense"] = tmrnadf.progress_apply(lambda row: check_UR_sense(row), axis=1)
    dict_ctx = create_dict_ctxs(tmrnadf)
    tmrnadf["tdna_name"] = tmrnadf.apply(lambda row: apply_nom(row, dict_ctx, nom_ext), axis=1)
    tmrnadf = modify_duplicate(tmrnadf)
    tmrnadf[["strain", "start", "end", "sense", "LCB_UP", "LCB_DOWN", "tdna_name"]].to_csv(
        f"{output_folder}/tdna_scheme_{nom_ext}.csv")
    create_genbanks(list_files_fasta, tmrnadf, file_ext)
    shutil.rmtree(f"{output_folder}/tmp/")
    
