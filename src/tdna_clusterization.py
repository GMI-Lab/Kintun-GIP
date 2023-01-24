#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 18:34:12 2020
for Synerclust
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

# Dictionary with anticodon codes
dct_ant = dict(ACC='gly4', ATG='his2', ACA='cys2', ACG='arg1', ATC='asp2',
               ATA='tyr2', AGG='pro4', CCT='arg2', AGC='ala3', AGA='ser6',
               ATT='asn2', CTG='gln2', CTC='glu2', CCG='arg4', AGT='thr4',
               CCA='trp1', CCC='gly1', TAT='ile2', GGT='thr1', CGA='ser3',
               CGC='ala4', CGG='pro1', GGG='pro2', TAG='leu3', GGA='ser2',
               TAA='leu2', GGC='ala1', TAC='val1', CGT='thr2', GTA='tyr1',
               GTC='asp1', GTG='his1', TTC='glu1', GTT='asn1', TTG='gln1',
               AAG='leu6', AAA='phe2', AAC='val4', TTT='lys1', CAT='met1',
               AAT='ile3', CAC='val3', CAA='leu5', CAG='leu4', TGT='thr3',
               TCT='arg3', GAT='ile1', CTT='lys2', TGC='ala2', TGA='ser1',
               TGG='pro3', GAG='leu1', TCG='arg5', GAC='val2', TCC='gly3',
               GAA='phe1', TCA='sec1', GCA='cys1', GCC='gly2', GCG='arg6',
               GCT='ser4', ACT='ser5', CTA='Sup1', TTA='Sup2', ssrA='ssrA')
# Amino Acids one letter code
aa_cod = {'tRNA-Cys': 'C', 'tRNA-Asp': 'D', 'tRNA-Ser': 'S',
          'tRNA-Gln': 'Q', 'tRNA-Lys': 'K', 'tRNA-Ile': 'I',
          'tRNA-Pro': 'P', 'tRNA-Thr': 'T', 'tRNA-Phe': 'F',
          'tRNA-Asn': 'N', 'tRNA-Gly': 'G', 'tRNA-His': 'H',
          'tRNA-Leu': 'L', 'tRNA-Arg': 'R', 'tRNA-Trp': 'W',
          'tRNA-Ala': 'A', 'tRNA-Val': 'V', 'tRNA-Glu': 'E',
          'tRNA-Tyr': 'Y', 'tRNA-Met': 'M', 'tRNA-Sec': 'U',
          'tRNA-Pyl': 'O'}

# Functions
def apply_ant_dict(ID_info):
    tmdna_id = ""
    if "Name=tmRNA" in ID_info:
        tmdna_id = "ssrA"
    else:
        anticodon = ID_info.split("(")[1].split(")")[0]
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
                
def write_fasta_and_clusterize(input_folder, output_folder, df,win1,win2):
    all_tdna_nom = {}
    for group_name, df_group in df.groupby(["nom_ant"]):
        # Create fasta
        fasta_file = output_folder + group_name + "_all.fasta"
        with open(fasta_file, "w") as salida:
            for row_index, row in df_group.iterrows():
                in_file = input_folder + row["strain"] + ".fasta"
                if row['sense'] == "+":
                    with open(in_file,"r") as entrada:
                        for record in SeqIO.parse(entrada,"fasta"):
                            DR_seq = str(record.seq[row["start"]-win1:row["start"]-win2]).upper()
                else:
                    with open(in_file,"r") as entrada:
                        for record in SeqIO.parse(entrada,"fasta"):
                            DR_seq = str(record.seq[row["end"]+win2:row["end"]+win1].reverse_complement()).upper()
                salida.write(">"+str(row["tdna_ind"])+"_-_"+row["nom_ant"]+"_-_"+row["strain"]+"\n")
                salida.write(DR_seq+"\n")
        
        # Clusterize with MMSeqs
        cmd = f"mmseqs easy-cluster {fasta_file} {fasta_file}_clusterRes_{str(win1)} {output_folder}/tmp/ -c 0.95 -v 0"
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

def create_genbanks(input_folder,output_folder,df):
    for group_name, df_group in df.groupby(["strain"]):
        with open(f"{input_folder}/{group_name}.fasta", "r") as record_file:
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
                    
                    # create new feature
                    new_feature = SeqFeature(
                        feature_location,
                        type=feature_type,
                        strand=feature_strand,
                        id=row["tdna_ind"],
                        qualifiers={
                           "locus_tag" : row["tdna_nom_cor"],
                           "product": row["ID"],
                           "Kintun_VLI-tDNAclass" : row["tdna_nom_cor"]
                           }
                        )
                    record.features.append(new_feature)
                    
                    for index,value in enumerate(row["corr_dists_dr"]):
                        if value[0] < value[1]:
                            dr_start_pos = ExactPosition(value[0])
                            dr_end_pos = ExactPosition(value[1])
                        else:
                            dr_start_pos = ExactPosition(value[1])
                            dr_end_pos = ExactPosition(value[0])
                        dr_feature_location = FeatureLocation(dr_start_pos, dr_end_pos)
                        dr_feature_type = "misc_feature"
                        id_string = f'DR_{row["tdna_nom_cor"]}_{index}'
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
                        
                    for index,value in enumerate(row["corr_dists_ur"]):
                        if value[0] < value[1]:
                            ur_start_pos = ExactPosition(value[0])
                            ur_end_pos = ExactPosition(value[1])
                        else:
                            ur_start_pos = ExactPosition(value[1])
                            ur_end_pos = ExactPosition(value[0])
                        ur_feature_location = FeatureLocation(ur_start_pos, ur_end_pos)
                        ur_feature_type = "misc_feature"
                        id_string = f'UR_{row["tdna_nom_cor"]}_{index}'
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
                with open(f"{output_folder}/{group_name}_VLIresults.gb", "w") as output_file:
                    SeqIO.write(record, output_file, "gb")

def check_exclusion(df,input_folder,output_folder):
    """
    A single chromosome cannot have
    two t(m)DNAS with the same name
    """
    repeated_tdnas = []
    corrected_tdna_nom = {}
    new_distances = []
    corrected_tdna_nom_distance = {}

    for group_name, df_group in df.groupby(["strain"]):
        if len(list(df_group["tdna_nom"])) != len(set(list(df_group["tdna_nom"]))):
            tdna_list = list(df_group["tdna_nom"])
            d = Counter(tdna_list)
            repeated_tdnas = repeated_tdnas + [item for item in d if d[item]>1]

    for tdna in list(dict.fromkeys(repeated_tdnas)):
        for group_name, df_group in df.groupby(["strain"]):
            sub_df = df_group.loc[df_group['nom_ant'] == tdna[:-1]]
            nom_dict = {}
            x = 2001
            while len(set(list(nom_dict.values()))) < len(sub_df):
                # deletes old file
                try:
                    filePath = f"{output_folder}/{tdna[:-1]}_all.fasta_clusterRes_{x-1000}_rep_seq.fasta"
                    os.remove(filePath)
                except FileNotFoundError:
                    pass
                # repeat UR analysis
                nom_dict = write_fasta_and_clusterize(input_folder,output_folder,sub_df,x,1)
                x += 2000
            try:
                filePath = f"{output_folder}/{tdna[:-1]}_all.fasta_clusterRes_{x-1000}_rep_seq.fasta"
                os.remove(filePath)
            except FileNotFoundError:
                pass
            new_distances.append((tdna,x))

    for tdna in list(dict.fromkeys(repeated_tdnas)):
        corrected_distance = np.max([i[1] for i in new_distances if i[0] == tdna])
        corrected_tdna_nom_distance[tdna] = corrected_distance
        sub_df = df.loc[df['nom_ant'] == tdna[:-1]]
        corrected_tdna_nom.update(write_fasta_and_clusterize(input_folder,output_folder,sub_df,corrected_distance,1))
    
    return corrected_tdna_nom,corrected_tdna_nom_distance


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

def conserved_downstream_blocks(list_files,input_folder,output_folder,df,win1):
    all_tdna_dists = {}
    for group_name, df_group in df.groupby(["tdna_nom_cor"]):
        if len(df_group) > 1:
            # ids in a list
            records_ids = []                        
            # Create fasta with DR regions
            fasta_file = output_folder + "dr_" + group_name + "_all.fasta"
            with open(fasta_file, "w") as salida:
                for row_index, row in df_group.iterrows():
                    strain_fasta = input_folder + row["strain"] + ".fasta"
                    if row['sense'] == "+":
                        DR_seq = extract_sequence_dr(strain_fasta, row["end"]+1, row["end"]+win1, "+")
                    else:
                        DR_seq = extract_sequence_dr(strain_fasta, row["start"]-win1, row["start"]-1, "-")
                    #salida
                    salida.write(">"+str(row["tdna_ind"])+"_-_"+row["tdna_nom_cor"]+"_-_"+row["strain"]+"\n")
                    salida.write(str(DR_seq)+"\n")
                    records_ids.append(row["tdna_ind"])
                                    
            # SibeliaZ and maf2synteny exec
            cmd_sibeliaz = f"sibeliaz -n -o {output_folder}/tmp_sibelia/ {fasta_file} ; maf2synteny -b 5000 -o {output_folder}/tmp_sibelia/ {output_folder}/tmp_sibelia/blocks_coords.gff"
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

def create_tdnas_scheme(input_folder, list_reps, df):
    with open("tdna_scheme.tsv","w") as output:
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
            with open(f"{input_folder}/{strain}.fasta","r") as strain_seq:
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
            with open("tdna_scheme.tsv","a") as output:
                output.write(f"{tdna_name}\t{prevalence}\t{strain}\t{ur_seq}\t{tdna_seq}\t{dr_seq}\n")

def tdna_clusterization(input_folder, output_folder, file_ext, nom_ext):
    # Create list of files .aragorn
    print(input_folder)
    list_files = glob.glob(f'{output_folder}/*/*.aragorn', recursive=True)
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
    nom_dict = write_fasta_and_clusterize(input_folder,output_folder,df,1001,1)
    df["tdna_nom"] = df.apply(lambda row : apply_new_nom(row["tdna_ind"],nom_dict), axis=1)

    # check exclusion nomenclature
    cor_nom_dict,cor_nom_dict_dist = check_exclusion(df,input_folder,output_folder)
    df["tdna_nom_cor"] = df.apply(lambda row : correct_nom(row["tdna_ind"],nom_dict,cor_nom_dict), axis=1)

    # UR sizes 
    df["uncorr_dists_ur"] = df.apply(lambda row : apply_dists_ur(row["nom_ant"],cor_nom_dict_dist), axis=1)
    df["corr_dists_ur"] = df.apply(lambda row : correct_distances_ur(row["sense"], row["start"], row["end"], row["uncorr_dists_ur"]), axis=1)

    #Calculate conserved downstream region
    dist_dr_cons = conserved_downstream_blocks(list_files,input_folder,output_folder,df,200000)
    df["uncorr_dists_dr"] = df.apply(lambda row : apply_dists(row["tdna_ind"],dist_dr_cons), axis=1)
    df["corr_dists_dr"] = df.apply(lambda row : correct_distances(row["sense"], row["start"], row["end"], row["uncorr_dists_dr"]), axis=1)

    #Export data to CSV file
    df.to_csv('export_tdnas_test.csv', index=False)

    #Create tDNAs scheme
    list_reps = glob.glob(f'{output_folder}/*_rep_seq.fasta', recursive=True)
    create_tdnas_scheme(input_folder,list_reps, df)

    #Create genbanks with annotations
    create_genbanks(input_folder,output_folder,df)
    
    #Remove cluster files
    for j in glob.glob(f'{output_folder}/*rep_seq.fasta', recursive=True):
        os.remove(j)
    shutil.rmtree(f"{output_folder}/tmp/")

