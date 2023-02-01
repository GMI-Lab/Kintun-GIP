#!/usr/bin/env python
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
import string
import subprocess
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import FeatureLocation
from collections import Counter
import numpy as np
import logging

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
               SERGCT='ser4', SERACT='ser5', SUPCTA='Sup1', SUPTTA='Sup2', 
               ILE2CAT='ile3', FMETCAT='fmet1', ssrA='ssrA')

# Functions
def apply_ant_dict(ID_info):
    tmdna_id = ""
    if "Name=tmRNA" in ID_info:
        tmdna_id = "ssrA"
    else:
        anticodon = ID_info.split("Name=")[1].replace("-","")
        tmdna_id = dct_ant[anticodon.upper()]
    return tmdna_id

def calculate_locus(infile):
    locus_id = 0
    list_tmdnas = []
    with open(infile, "r") as entrada:
        for line in entrada.readlines():
            sep = line.strip().split("\t")
            list_tmdnas.append((sep[3],sep[4],sep[6]))
    return list_tmdnas
                
def write_fasta_and_clusterize(input_folder, file_ext, output_folder, df, win1 ,win2, threads):
    all_tdna_nom = {}
    for group_name, df_group in df.groupby(["nom_ant"]):
        # Create fasta
        fasta_file = f"{output_folder}/{group_name}_all.fasta"
        with open(fasta_file, "w") as salida:
            for row_index, row in df_group.iterrows():
                in_file = f"{input_folder}/{row['strain']}.{file_ext}"
                if row['sense'] == "+":
                    with open(in_file,"r") as entrada:
                        for record in SeqIO.parse(entrada,"fasta"):
                            DR_seq = str(record.seq[row["start"]-win1:row["start"]-win2]).upper()
                else:
                    with open(in_file,"r") as entrada:
                        for record in SeqIO.parse(entrada,"fasta"):
                            DR_seq = str(record.seq[row["end"]+win2:row["end"]+win1].reverse_complement()).upper()
                salida.write(f">{str(row['tdna_ind'])}_-_{row['nom_ant']}_-_{row['strain']}\n")
                salida.write(f"{DR_seq}\n")
        # Clusterize with MMSeqs
        cmd = f"mmseqs easy-cluster {fasta_file} {fasta_file}_clusterRes_{str(win1)} {output_folder}/tmp/ " \
              f"-c 0.95 -v 0 --threads {threads}"
        p = subprocess.run(cmd, shell=True)
        #os.remove(fasta_file)
        os.remove(f"{fasta_file}_clusterRes_{str(win1)}_all_seqs.fasta")
        #Abre el archivo de las asignaciones
        with open(f"{fasta_file}_clusterRes_{str(win1)}_cluster.tsv", "r") as tsv_file:
            #Lee el archivo de asignaciones
            lines = [i.strip().split("\t") for i in tsv_file.readlines()]
            #Procesa los pares asociados
            pairs = []
            for i in lines:
                pairs.append((i[0],i[1]))
            result = {}
            for sub_tuple in pairs:
               # checking the first element of the tuple in the result
               if sub_tuple[0] in result:
                  # adding the current tuple values without first one
                  result[sub_tuple[0]] = (*result[sub_tuple[0]], *sub_tuple[1:])
               else:
                  # adding the tuple
                  result[sub_tuple[0]] = sub_tuple
            # en groups estan los clusters
            groups = [tuple(set(i)) for i in result.values()]
            groups.sort(key=lambda t: len(t), reverse=True)
            alphabet = list(string.ascii_uppercase) + [
                i + j for j in list(string.ascii_uppercase) for i in list(string.ascii_uppercase)]
            tdna_nom = {}
            for index, value in enumerate(groups):
                for tdna in value:
                    tdna_nom[tdna.split("_-_")[0]] = tdna.split("_-_")[1] + alphabet[index]
            all_tdna_nom.update(tdna_nom)
            os.remove(f"{fasta_file}_clusterRes_{str(win1)}_cluster.tsv")
            os.remove(fasta_file)
    return all_tdna_nom

def apply_new_nom(tdna_id, dict_nom):
    return dict_nom[str(tdna_id)]

def correct_nom(tdna_id, dict_nom, cor_dict_nom):
    try:
        return cor_dict_nom[str(tdna_id)]
    except KeyError:
        return dict_nom[str(tdna_id)]

def create_genbanks(input_folder,file_ext,output_folder,df,prefix):
    for group_name, df_group in df.groupby(["strain"]):
        with open(f"{input_folder}/{group_name}.{file_ext}", "r") as record_file:
            for record in SeqIO.parse(record_file,"fasta"):
                for row_index, row in df_group.iterrows():
                    # positions
                    start_pos = ExactPosition(row["start"])
                    end_pos = ExactPosition(row["end"])
                    feature_location = FeatureLocation(start_pos, end_pos)
                    # type
                    if "Name=tmRNA" in row["ID"]:
                        feature_type = "tmRNA"
                    else:
                        feature_type = "tRNA"
                    # strand
                    if row["sense"] == "+":
                        feature_strand = 1
                    else:
                        feature_strand = -1
                    # create new tmrna feature
                    new_feature = SeqFeature(
                        feature_location,
                        type=feature_type,
                        strand=feature_strand,
                        id=row["tdna_ind"],
                        qualifiers={
                           "locus_tag" : f"{prefix}_{row['tdna_nom_cor']}",
                           "product": row["ID"],
                           "Kintun_GPI-tmDNAclass" : f"{prefix}_{row['tdna_nom_cor']}"
                           }
                        )
                    record.features.append(new_feature)
                    # create new dr feature
                    for index,value in enumerate(row["corr_dists_dr"]):
                        if value[0] < value[1]:
                            dr_start_pos = ExactPosition(value[0])
                            dr_end_pos = ExactPosition(value[1])
                        else:
                            dr_start_pos = ExactPosition(value[1])
                            dr_end_pos = ExactPosition(value[0])
                        dr_feature_location = FeatureLocation(dr_start_pos, dr_end_pos)
                        dr_feature_type = "misc_feature"
                        id_string = f'DR_{prefix}_{row["tdna_nom_cor"]}_{index}'
                        dr_feature = SeqFeature(
                        dr_feature_location,
                        type = dr_feature_type,
                        strand = 0,
                        id = id_string,
                        qualifiers={
                           "locus_tag" : id_string,
                           }
                        )
                        record.features.append(dr_feature)
                    # create new ur feature
                    for index,value in enumerate(row["corr_dists_ur"]):
                        if value[0] < value[1]:
                            ur_start_pos = ExactPosition(value[0])
                            ur_end_pos = ExactPosition(value[1])
                        else:
                            ur_start_pos = ExactPosition(value[1])
                            ur_end_pos = ExactPosition(value[0])
                        ur_feature_location = FeatureLocation(ur_start_pos, ur_end_pos)
                        ur_feature_type = "misc_feature"
                        id_string = f'UR_{prefix}_{row["tdna_nom_cor"]}_{index}'
                        ur_feature = SeqFeature(
                        ur_feature_location,
                        type = ur_feature_type,
                        strand = 0,
                        id = id_string,
                        qualifiers={
                           "locus_tag" : id_string,
                           }
                        )
                        record.features.append(ur_feature)      
                    record.annotations["molecule_type"] = "DNA"
                    record.annotations["form"] = "double stranded"
                    record.annotations["topology"] = "circular"
                # export genbank
                with open(f"{output_folder}/{group_name}_kintun-clust.gb", "w") as output_file:
                    SeqIO.write(record, output_file, "gb")

def detect_repeated_names(df,tdna,x):
    check_list = []
    for group_name, df_group in df.groupby(["strain"]):
        if x == 1001:
            if len(list(df_group[f"tdna_nom"])) != len(set(list(df_group[f"tdna_nom"]))):
                check_list.append(group_name)
            else:
                continue
        else:
# tdna_nom_ile1_2001
            if len(list(df_group[f"tdna_nom_{tdna}_{str(x)}"])) != len(set(list(df_group[f"tdna_nom_{tdna}_{str(x)}"]))):
                check_list.append(group_name)
    return check_list

def check_exclusion(df,input_folder,file_ext,output_folder,threads):
    """
    A single chromosome cannot have
    two t(m)DNAS with the same name
    """
    repeated_tdnas_prev = []
    new_distances = {}
    corrected_tdna_nom_distance = {}
    corrected_tdna_nom = {}

    for group_name, df_group in df.groupby(["strain"]):
        if len(list(df_group["tdna_nom"])) != len(set(list(df_group["tdna_nom"]))):
            tdna_list = list(df_group["tdna_nom"])
            d = Counter(tdna_list)
            repeated_tdnas_prev = repeated_tdnas_prev + [item for item in d if d[item]>1]
    repeated_tdnas = []
    for i in repeated_tdnas_prev:
        if i[:-1] not in repeated_tdnas:
            repeated_tdnas.append(i[:-1])
        else:
            continue
    logging.info(f"Detected problems with the following tmDNAs {' '.join(repeated_tdnas)}")
    
    for tdna in repeated_tdnas:
        df_tdna = pd.DataFrame()
        for group_name, df_group in df.groupby(["nom_ant"]):
            if group_name == tdna:
                df_tdna = pd.concat([df_tdna,df_group])
        df_tdna.to_csv(f'{output_folder}/export_{tdna}_tdnas_test_prev.csv', index=False)
        # First try
        nom_dict_new = write_fasta_and_clusterize(input_folder, file_ext, output_folder, df_tdna, 251, 1, threads)
        logging.info(f"Creating tdna_nom_{tdna}_{251} column")
        df_tdna[f"tdna_nom_{tdna}_{251}"] = df_tdna.apply(lambda row: apply_new_nom(row["tdna_ind"], nom_dict_new), axis=1)
        df_tdna.to_csv(f'{output_folder}/export_{tdna}_tdnas_test_{251}.csv', index=False)
        # Loop
        window_values = [251,500,2001,4001,6001,8001,10001,20001,50001,100001]
        for index,value in enumerate(window_values):
            if index == 0:
                problem_detected = detect_repeated_names(df_tdna,tdna,value)
                logging.info(f"index = {index}, found problem with {' '.join(problem_detected)}")
            else:
                problem_detected = detect_repeated_names(df_tdna,tdna,window_values[index-1])
                logging.info(f"index = {index}, found problem with {' '.join(problem_detected)}")
            if len(problem_detected) > 1:
                logging.info(f"Checking {tdna} with {str(value)} pb windows")
                nom_dict_new = write_fasta_and_clusterize(input_folder, file_ext, output_folder, df_tdna, value, 1, threads)
                logging.info(f"Creating tdna_nom_{tdna}_{str(value)} column")
                df_tdna[f"tdna_nom_{tdna}_{str(value)}"] = df_tdna.apply(lambda row: apply_new_nom(row["tdna_ind"], nom_dict_new), axis=1)
                df_tdna.to_csv(f'{output_folder}/export_{tdna}_tdnas_test_{str(value)}.csv', index=False)
                if len(detect_repeated_names(df_tdna,tdna,value)) > 1:
                    continue
                else:
                    new_distances[tdna] = value
                    break
            else:
                new_distances[tdna] = value
                break
        try:
            new_distances[tdna]
        except KeyError:
            new_distances[tdna] = 1001       

    for tdna in repeated_tdnas:
        corrected_distance = new_distances[tdna]
        corrected_tdna_nom_distance[tdna] = corrected_distance
        sub_df = df.loc[df['nom_ant'] == tdna]
        corrected_tdna_nom.update(write_fasta_and_clusterize(input_folder, file_ext, output_folder, sub_df, corrected_distance, 1, threads))
    return corrected_tdna_nom, corrected_tdna_nom_distance

def extract_sequence_dr(fasta_file, start, end, sense):
    with open(fasta_file, "r") as entrada:
        for record in SeqIO.parse(entrada,"fasta"):
            # concatenates_chromosomes
            seq_sum = record.seq + record.seq
            # new_coordinates
            if start < 0 and end < 0:
                dr_start = len(record.seq) + start
                dr_end = len(record.seq) + end
            elif start < 0 and end >= 0:
                dr_start = len(record.seq) + start
                dr_end = len(record.seq) + end
            else:
                dr_start = start
                dr_end = end
            # chech_sense and export
            if sense == "+":
                dr_seq = seq_sum[dr_start:dr_end].upper()
            else:
                dr_seq = seq_sum[dr_start:dr_end].reverse_complement()
            # return
            return dr_seq

def conserved_downstream_blocks(list_files,input_folder,file_ext,output_folder,df,win1,threads):
    all_tdna_dists = {}
    for group_name, df_group in df.groupby(["tdna_nom_cor"]):
        if len(df_group) > 1:
            # ids in a list
            records_ids = []                        
            # Create fasta with DR regions
            fasta_file = output_folder + "dr_" + group_name + "_all.fasta"
            with open(fasta_file, "w") as salida:
                for row_index, row in df_group.iterrows():
                    strain_fasta = f"{input_folder}/{row['strain']}.{file_ext}"
                    if row['sense'] == "+":
                        DR_seq = extract_sequence_dr(strain_fasta, row["end"]+1, row["end"]+win1, "+")
                    else:
                        DR_seq = extract_sequence_dr(strain_fasta, row["start"]-win1, row["start"]-1, "-")
                    #salida
                    salida.write(">"+str(row["tdna_ind"])+"_-_"+row["tdna_nom_cor"]+"_-_"+row["strain"]+"\n")
                    salida.write(str(DR_seq)+"\n")
                    records_ids.append(row["tdna_ind"])                                    
            # SibeliaZ and maf2synteny exec
            cmd_sibeliaz = f"sibeliaz -n -t {threads} -o {output_folder}/tmp_sibelia/ {fasta_file} ; maf2synteny -b 5000 -o {output_folder}/tmp_sibelia/ {output_folder}/tmp_sibelia/blocks_coords.gff"
            p = subprocess.run(cmd_sibeliaz, shell=True)               
            with open(f"{output_folder}/tmp_sibelia/5000/blocks_coords.txt","r") as lcbs_file:
                fields = "".join(lcbs_file.readlines()).split(f"{80*'-'}\n")
                # Extract information
                seqs = [ i.strip().split("\n") for i in fields if not "Block #" in i]
                groups = [ i.strip().split("\n") for i in fields if "Block #" in i]
                conserved_lcb = [i for i in groups if len(i) == len(list_files)+2]
                # Obtain coords
                dct_codes = {}
                lcbs = []
                for i in seqs[0][1:]:
                    line = i.strip().split("\t")
                    dct_codes[line[0]] = line[-1]
                for block in conserved_lcb:
                    for rline in block[2:]:
                        lcb = rline.strip().split("\t")
                        lcbs.append((dct_codes[lcb[0]],int(lcb[-3]),int(lcb[-2])))
                            
            # Process data
            tdna_nom = {}
            for index, value in enumerate(records_ids):
                tdna_nom[value] = [(i[1],i[2]) for i in lcbs if i[0].split("_-_")[0] == str(value)]
                all_tdna_dists.update(tdna_nom)
            shutil.rmtree(f"{output_folder}/tmp_sibelia/")
            os.remove(fasta_file)
        elif len(df_group) == 1:
            for row_index, row in df_group.iterrows():
                all_tdna_dists[row["tdna_ind"]] = [(1,2001)]
            #eliminar la carpeta
    return all_tdna_dists

def apply_dists(tdna_id, dict_dist):
    return dict_dist[tdna_id]

def apply_dists_ur(nom_ant, dict_dist):
    dict_dist_ant = {}
    distances = [dict_dist[i] for i in dict_dist.keys()]
    nom_ant_cor = [i[:-1] for i in dict_dist.keys()]
    for index,value in enumerate(nom_ant_cor):
        dict_dist_ant[value] = distances[index]
    if nom_ant in dict_dist_ant.keys():
        return dict_dist_ant[nom_ant]
    else:
        return 1001

def correct_distances(sense, tdna_start, tdna_end, coordinates):
    corrected_coordinates = []
    if sense == "+":
        for i in coordinates:
            if i[1] < 2001:
                pos_start = i[0]
                pos_end = i[1]
                new_start = tdna_end + pos_start
                new_end = tdna_end + pos_end
                corrected_coordinates.append((new_start,new_end))
            elif i[1] > 2001:
                pos_start = i[0]
                pos_end = i[0] + 2001 
                new_start = tdna_end + pos_start
                new_end = tdna_end + pos_end
                corrected_coordinates.append((new_start,new_end))
            else:
                corrected_coordinates.append((tdna_end+1,tdna_end+2001))
    else:
        for i in coordinates:
            if i[1] < 2001:
                pos_start = i[0]
                pos_end = i[1]
                new_start = tdna_start - pos_start
                new_end = tdna_start - pos_end
                corrected_coordinates.append((new_start,new_end))
            elif i[1] > 2001:
                pos_start = i[0]
                pos_end = i[0] + 2001
                new_start = tdna_start - pos_start
                new_end = tdna_start - pos_end
                corrected_coordinates.append((new_start,new_end))
                break
            else:
                corrected_coordinates.append((tdna_start-1,tdna_start-2001))
    return corrected_coordinates

def correct_distances_ur(sense, tdna_start, tdna_end, distance):
    corrected_coordinates = []
    if sense == "-":
        new_start = tdna_end + 1
        new_end = tdna_end + distance
        corrected_coordinates.append((new_start,new_end))
    else:
        new_start = tdna_start - distance
        new_end = tdna_start - 1
        corrected_coordinates.append((new_start,new_end))
    return corrected_coordinates

def create_tdnas_scheme(input_folder,file_ext,output_folder,list_reps, df, prefix):
    with open(f"{output_folder}/{prefix}_tdna_scheme.tsv","w") as output:
        output.write("#tdna_name\tprevalence\tstrain\tur_seq\ttdna_seq\tdr_seq\n")            
    # check by anticodon
    for file in list_reps:
        ids = []
        with open(file,"r") as in_file:
            for line in in_file:
                if line[0] != ">":
                    continue
                else:
                    ids.append(line.strip().replace(">","").split("_-_"))        
        for tdna_class in ids:
            rep_locus = df.loc[df["tdna_ind"] == int(tdna_class[0])]
            tdna_name = rep_locus.tdna_nom_cor.item()
            strain = rep_locus.strain.item()
            prevalence = len(df.loc[df.tdna_nom_cor == tdna_name])
            with open(f"{input_folder}/{strain}.{file_ext}","r") as strain_seq:
                for record in SeqIO.parse(strain_seq,"fasta"):
                    if rep_locus.sense.item() == "+":
                        tdna_seq = record.seq[rep_locus.start.item():rep_locus.end.item()+1]
                        ur_seq = ''
                        for j in rep_locus.corr_dists_ur.item():
                            ur_seq += record.seq[j[0]:j[1]+1]
                        dr_seq = ''
                        for k in rep_locus.corr_dists_dr.item():
                            dr_seq += record.seq[k[0]:k[1]+1]
                    else:
                        tdna_seq = record.seq[rep_locus.start.item():rep_locus.end.item()+1].reverse_complement()
                        ur_seq = ''
                        for j in rep_locus.corr_dists_ur.item():
                            ur_seq += record.seq[j[0]:j[1]+1].reverse_complement()
                        dr_seq = ''
                        for k in rep_locus.corr_dists_dr.item():
                            dr_seq += record.seq[k[1]:k[0]+1].reverse_complement()
            with open(f"{output_folder}/{prefix}_tdna_scheme.tsv","a") as output:
                output.write(f"{prefix}_{tdna_name}\t{prevalence}\t{strain}\t{ur_seq}\t{tdna_seq}\t{dr_seq}\n")

def tdna_clusterization(input_folder,output_folder,file_ext,nom_ext,threads):
    # Create list of files .aragorn
    list_files = glob.glob(f'{output_folder}/*/*.trnascanse', recursive=True)
    # Create empty dataframes with column numbers
    df = pd.DataFrame(columns=["strain","method","type","start","end","na1","sense","na2","ID"])
    # For each file add data to the df Dataframe
    for infile in list_files:
        df2 = pd.read_csv(infile, sep="\t", names=["strain","method","type","start","end","na1","sense","na2","ID"])
        df = pd.concat([df,df2])
    #Delete no data columns
    df.drop(["method","type","na1","na2"], axis=1, inplace=True)
    #Corrects indexes
    df.reset_index(inplace=True, drop=True)
    #Add some info
    df["nom_ant"] = df.apply(lambda row : apply_ant_dict(row["ID"]), axis=1)
    df['tdna_ind'] = range(1, len(df) + 1)
    # defines nomenclature
    nom_dict = write_fasta_and_clusterize(input_folder,file_ext,output_folder,df,1001,1,threads)
    df["tdna_nom"] = df.apply(lambda row : apply_new_nom(row["tdna_ind"],nom_dict), axis=1)
    # check exclusion nomenclature
    cor_nom_dict,cor_nom_dict_dist = check_exclusion(df,input_folder,file_ext,output_folder,threads)
    df["tdna_nom_cor"] = df.apply(lambda row : correct_nom(row["tdna_ind"],nom_dict,cor_nom_dict), axis=1)
    # UR sizes 
    df["uncorr_dists_ur"] = df.apply(lambda row : apply_dists_ur(row["nom_ant"],cor_nom_dict_dist), axis=1)
    df["corr_dists_ur"] = df.apply(lambda row : correct_distances_ur(row["sense"], row["start"], row["end"], row["uncorr_dists_ur"]), axis=1)
    #Calculate conserved downstream region
    dist_dr_cons = conserved_downstream_blocks(list_files,input_folder,file_ext,output_folder,df,200000,threads)
    df["uncorr_dists_dr"] = df.apply(lambda row : apply_dists(row["tdna_ind"],dist_dr_cons), axis=1)
    df["corr_dists_dr"] = df.apply(lambda row : correct_distances(row["sense"], row["start"], row["end"], row["uncorr_dists_dr"]), axis=1)
    #Create tDNAs scheme
    list_reps = glob.glob(f'{output_folder}/*_rep_seq.fasta', recursive=True)
    create_tdnas_scheme(input_folder,file_ext,output_folder,list_reps, df, nom_ext)
    #Create genbanks with annotations
    create_genbanks(input_folder,file_ext,output_folder,df,nom_ext)
    #Remove cluster files
    for j in glob.glob(f'{output_folder}/*rep_seq.fasta', recursive=True):
        os.remove(j)
    #Remove tmp folder
    shutil.rmtree(f"{output_folder}/tmp/")

   
