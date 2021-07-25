#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kintun-VLI
by Pablo Lorca-Orloff, Carlos Serrano and Camilo Berrios
August 2020
v1.0
"""

from Bio import SeqIO
from scipy.cluster.hierarchy import fclusterdata
import numpy as np
import string
import multiprocessing
from functools import partial

# DEFINE DICTIONARIES FOR NOMENCLATURE BASIS

def create_dics():
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
    aa_cod = {'tRNA-Cys': 'C', 'tRNA-Asp': 'D', 'tRNA-Ser': 'S',
              'tRNA-Gln': 'Q', 'tRNA-Lys': 'K', 'tRNA-Ile': 'I',
              'tRNA-Pro': 'P', 'tRNA-Thr': 'T', 'tRNA-Phe': 'F',
              'tRNA-Asn': 'N', 'tRNA-Gly': 'G', 'tRNA-His': 'H',
              'tRNA-Leu': 'L', 'tRNA-Arg': 'R', 'tRNA-Trp': 'W',
              'tRNA-Ala': 'A', 'tRNA-Val': 'V', 'tRNA-Glu': 'E',
              'tRNA-Tyr': 'Y', 'tRNA-Met': 'M', 'tRNA-Sec': 'U',
              'tRNA-Pyl': 'O'}
    return dct_ant, aa_cod

# PART I. COLECT ALL t(m)DNAs loci
# Class locus have all information for each locus to be considerated
# to classification


class locus:
    def __init__(self, strain, tag):
        self.strain = strain
        self.tag = "%s_%s" % (strain, tag)
        self.sense = ''
        self.tdnas = []
        self.anticodon_in_locus = []
        self.UR_posI = ''
        self.UR_posII = ''
        self.UR_posIII = ''
        self.UR_groups = []
        self.pos_senses = []
        self.tdnas_in_UR = ''
        self.rdnas_in_UR = ''
        self.UR_I_class = []
        self.UR_II_class = []
        self.UR_III_class = []
        self.classification = []
        self.context_list = []
        self.record = []


def identifica_locus(strain_codes):
    global next_feature_index
    dct_ant, aa_cod = create_dics()
    loci_collection = []
    tag = 1
    for strain in strain_codes:
        # Open GenBank file with annotated sequences
        strain_sequence = SeqIO.read(strain, 'genbank')
        annotations_list = strain_sequence.features
        # Busca en las listas los tDNAs o tmDNAs (sí, también tmDNAs, uhlalá)
        while feature_index in range(0,len(annotations_list)):
            feature = annotations_list[feature_index]
            if feature.type in ['tRNA', 'tmRNA']:
                new_locus = locus(strain, tag)
                new_locus.sense = feature.strand
                new_locus.tdnas.append(feature)
                new_locus.record = feature
                for next_feature_index in range(feature_index+1, len(annotations_list)):
                    if annotations_list[next_feature_index].type == 'tRNA':
                        new_locus.tdnas.append(annotations_list[next_feature_index])
                        continue
                    else:
                        break
                locus_end = next_feature_index
                # Identifica si es un operon o una unidad mono o policistrónica
                if len([new_locus.tdnas]) > 1:
                    new_locus.classification = 'tdna_cluster'
                else:
                    new_locus.classification = 'tdna_single'
                # Identifica las tres posiciones río arriba
                features_ur, tdnas_ur, rdnas_ur = ([],) * 3
                features_dr, tdnas_dr, rdnas_dr = ([],) * 3
                for feature in (annotations_list[:feature_index][::-1] + annotations_list[::-1]):
                    if len(features_ur) < 3:
                        if feature.type == 'CDS':
                            features_ur.append(feature)
                        elif feature.type == 'tRNA' and 'tRNA-Stop' not in feature.qualifiers['note'][0]:
                            tdnas_ur.append(
                                aa_cod[feature.qualifiers['product'][0]])
                        elif feature.type == 'rRNA':
                            rdnas_ur.append(feature.qualifiers['product'][0])
                        else:
                            continue
                    else:
                        break
                for feature in annotations_list[next_feature_index:] + annotations_list:
                    if len(features_dr) < 3:
                        if feature.type == 'CDS':
                            features_dr.append(feature)
                        elif feature.type == 'tRNA' and 'tRNA-Stop' not in feature.qualifiers['note'][0]:
                            tdnas_dr.append(
                                aa_cod[feature.qualifiers['product'][0]])
                        elif feature.type == 'rRNA':
                            rdnas_dr.append(feature.qualifiers['product'][0])
                        else:
                            continue
                    else:
                        break
                position_keys = ['I', 'II', 'III']
                for k in [0, 1, 2]:
                    if new_locus.sense == 1:
                        # UR features, tdnas and rdnas
                        exec(
                            "new_locus.UR_pos%s = features_ur[%s]" %
                            (position_keys[k], str(k))
                            )
                        new_locus.tdnas_in_UR = tdnas_ur
                        new_locus.rdnas_in_UR = rdnas_ur
                        # DR features, tdnas, rdnas
                        exec(
                            "new_locus.DR_pos%s = features_dr[%s]" %
                            (position_keys[k], str(k))
                            )
                        new_locus.tdnas_in_DR = tdnas_dr
                        new_locus.rdnas_in_DR = rdnas_dr
                    else:
                        # UR features, tdnas and rdnas
                        exec(
                            "new_locus.UR_pos%s = features_dr[%s]" %
                            (position_keys[k], str(k))
                        )
                        new_locus.tdnas_in_UR = tdnas_dr
                        new_locus.rdnas_in_UR = rdnas_dr
                        # DR features, tdnas, rdnas
                        exec(
                            "new_locus.DR_pos%s = features_ur[%s]" %
                            (position_keys[k], str(k))
                        )
                        new_locus.tdnas_in_DR = tdnas_ur
                        new_locus.rdnas_in_DR = rdnas_ur
                # Consider CDS sense in loci
                if new_locus.sense == 1:
                    for record in [
                            new_locus.UR_posI, new_locus.UR_posII,
                            new_locus.UR_posIII]:
                        if record.strand == 1:
                            sense_code = "_1"
                        else:
                            sense_code = "_0"
                        new_locus.UR_groups.append(
                            record.qualifiers['synerclust_group'][0] +
                            sense_code)
                    for record in [
                            new_locus.DR_posI, new_locus.DR_posII,
                            new_locus.DR_posIII]:
                        if record.strand == 1:
                            sense_code = "_1"
                        else:
                            sense_code = "_0"
                        new_locus.DR_groups.append(
                            record.qualifiers['synerclust_group'][0] +
                            sense_code)
                else:
                    for record in [
                            new_locus.UR_posI, new_locus.UR_posII,
                            new_locus.UR_posIII]:
                        if record.strand == 1:
                            sense_code = "_0"
                        else:
                            sense_code = "_1"
                        new_locus.UR_groups.append(
                            record.qualifiers['synerclust_group'][0] +
                            sense_code)
                    for record in [
                            new_locus.DR_posI, new_locus.DR_posII,
                            new_locus.DR_posIII]:
                        if record.strand == 1:
                            sense_code = "_0"
                        else:
                            sense_code = "_1"
                        new_locus.DR_groups.append(
                            record.qualifiers['synerclust_group'][0] +
                            sense_code)
                anticodons = []
                for element in new_locus.tdnas:
                    if element.type == 'tRNA':
                        anticodons.append(
                            element.qualifiers['note'][0].split(
                                '(')[1].strip(')').upper())
                    else:
                        anticodons.append('ssrA')
                new_locus.anticodon_in_locus = anticodons
                loci_collection.append(new_locus)
                feature_index = locus_end
                tag += 1
            else:
                feature_index += 1
                continue
    return loci_collection


def set_anticodons(loci_collection):
    ant_out = list(set([tdna.anticodon_in_locus for tdna in loci_collection]))
    return ant_out


def generate_binary_data_UR(anticodon, loci_to_analyze):
    # Collect CDS identifiers from loci to analyze
    posI_class =   [loci.UR_groups[0] for loci in loci_to_analyze]
    posII_class =  [loci.UR_groups[1] for loci in loci_to_analyze]
    posIII_class = [loci.UR_groups[2] for loci in loci_to_analyze]
    # list_categories have all the identifiers
    list_categories = list(set(posI_class + posII_class + posIII_class))
    # Create binary matrix for clustering
    binary_pos_collection = []
    for loci in loci_to_analyze:
        binary_loci = [1 if cat in loci.UR_groups else 0 for cat in list_categories]
        binary_pos_collection.append(binary_loci)
    return binary_pos_collection

def generate_binary_data_DR(anticodon, loci_to_analyze):
    # Collect CDS identifiers from loci to analyze
    posI_class = [loci.DR_groups[0] for loci in loci_to_analyze]
    posII_class = [loci.DR_groups[1] for loci in loci_to_analyze]
    posIII_class = [loci.DR_groups[2] for loci in loci_to_analyze]
    # list_categories have all the identifiers
    list_categories = list(set(posI_class + posII_class + posIII_class))
    # Create binary matrix for clustering
    binary_pos_collection = []
    for loci in loci_to_analyze:
        binary_loci = [1 if cat in loci.UR_groups else 0 for cat in list_categories]
        binary_pos_collection.append(binary_loci)
    return binary_pos_collection

def clustering_loci(ant_to_analyze, all_loci):
    dct_ant, aa_cod = create_dics()
    # Select loci
    for ant_cod in ant_to_analyze:
        coll_loci = [loci for loci in all_loci if ant_cod in loci.anticodon_in_locus]
        if len(coll_loci) <= 1:
            coll_loci[0].context_list.append("%s_%s_%s" % (dct_ant[ant_cod], '1', '1'))
        else:
            sequences_by_strain = {}
            for loci in coll_loci:
                sequences_by_strain[loci.strain] = []
            binary_data_ur = generate_binary_data_UR(ant_cod, coll_loci)
            labels_ur= fclusterdata(binary_data_ur, len(coll_loci),
                                    criterion='maxclust', metric='hamming',
                                    method='ward')
            for loci in coll_loci:
                sequences_by_strain[loci.strain] = []
            binary_data_dr = generate_binary_data_UR(ant_cod, coll_loci)
            labels_dr= fclusterdata(binary_data_dr, len(coll_loci),
                                    criterion='maxclust', metric='hamming',
                                    method='ward')
            for index in range(len(coll_loci)):
                coll_loci[index].context_list.append("%s_%s_%s" % (
                    dct_ant[ant_cod], labels_ur[index], labels_dr[index]))


class context_type:
    def __init__(self, code, num_strains, max_tdnas):
        self.code = code
        self.num_strains = num_strains
        self.max_tdnas = max_tdnas
        self.letters = []

    def __repr__(self):
        return repr((self.code, self.num_strains,
                     self.max_tdnas, self.letters))


def create_letter_list():
    list_alpha_original = list(string.ascii_uppercase)
    list_alpha = []
    for i in range(len(list_alpha_original)):
        list_alpha.append(list_alpha_original[i])
    for i in range(len(list_alpha_original)):
        for j in range(len(list_alpha_original)):
            list_alpha.append(list_alpha_original[i]+list_alpha_original[j])
    for i in range(len(list_alpha_original)):
        for j in range(len(list_alpha_original)):
            for k in range(len(list_alpha_original)):
                new_letter = (list_alpha_original[i] + list_alpha_original[j]
                              + list_alpha_original[k])
                list_alpha.append(new_letter)
    return list_alpha


def determine_nomenclature(all_loci, ant_to_analyze, prefix):
    dct_ant, aa_cod = create_dics()
    for ant_cod in ant_to_analyze:
        coll_loci = [loci for loci in all_loci if ant_cod in loci.anticodon_in_locus]
        contexts_types = []
        for locus in coll_loci:
            for context in locus.context_list:
                if dct_ant[ant_cod] in context and context not in contexts_types:
                    contexts_types.append(context)
                else:
                    continue
        contexts_types_info = []
        for context_type in contexts_types:
            strains = []
            list_number_tdnas = []
            for locus in coll_loci:
                if context in locus.context_list:
                    strains.append(locus.strain)
                    number_tdnas = []
                    for element in locus.tdnas:
                        if element.type == "tRNA" \
                                and element.qualifiers['note'][0].split('(')[1].strip(')').upper() == ant_cod:
                            number_tdnas.append(element)
                        elif element.type == "tmRNA":
                            number_tdnas.append(element)
                        else:
                            continue
                list_number_tdnas.append(len(number_tdnas))
            contexts_type = context_type(context, len(strains), np.max(list_number_tdnas))
            contexts_types_info.append(contexts_type)
        contexts_types_info = sorted(contexts_types_info, key=lambda context_type: context_type.max_tdnas, reverse=True)
        contexts_types_info = sorted(contexts_types_info, key=lambda context_type: context_type.num_strains,
                                     reverse=True)
        list_alpha = create_letter_list()
        #Cambié el sentido de esta lista, atencion
        context.letters = [list_alpha[i] for i in range(0,context.max_tdnas)]

        for locus in coll_loci:
            for loci_type in locus.context_list:
                if dct_ant[ant_cod] in loci_type:
                    locus_type = loci_type
                    break
            index_type = [context.code for context in contexts_types_info
                          ].index(loci_type)
            nomenclature_letters = contexts_types_info[index_type].letters
            tdnas_to_name = []
            for element in locus.tdnas:
                if element.type == "tRNA" and element.qualifiers['note'][0].split('(')[1].strip(')').upper() == ant_cod:
                    tdnas_to_name.append(element)
                elif element.type == "tmRNA":
                    tdnas_to_name.append(element)
            if locus.sense == 1:
                x = 0
                for tdna in tdnas_to_name:
                    tdna.qualifiers['Kintun_VLI-tDNAclass'] = [prefix + "_" + 
                                                     loci_type[:4] +
                                                     nomenclature_letters[x]]
                    x = x + 1
            else:
                x = 0
                for tdna in tdnas_to_name[::-1]:
                    tdna.qualifiers['Kintun_VLI-tDNAclass'] = [prefix + "_" +
                                                               loci_type[:4] +
                                                               nomenclature_letters[x]]
                    x = x + 1


def create_annotated_gbks(loci_collection, strain):
    new_features = []
    for locus in loci_collection:
        if locus.strain == strain:
            new_feature = [tdna for tdna in locus.tdnas]
            for tdna in new_feature:
                new_features.append(tdna)
        else:
            continue
    qualifiers_coll = []
    for feature in new_features:
        tdna_name = feature.qualifiers['Kintun_VLI-tDNAclass'][0]
        if tdna_name not in qualifiers_coll:
            qualifiers_coll.append(tdna_name)
        else:
            qualifiers_coll.append(tdna_name)
            feature.qualifiers['Kintun_VLI-tDNAclass'][0] = tdna_name + "_" + str(qualifiers_coll.count(tdna_name))
    with open(strain) as infile:
        seq_annot = SeqIO.read(infile, 'genbank')
        for feature in seq_annot.features:
            if feature.type not in ["tRNA", "tmRNA"]:
                new_features.append(feature)
            else:
                continue
    seq_annot.features = new_features
    with open(strain.replace("_SC_results.gbk",
                             "_KintunVLI_results.gbk"), 'w') as salida:
        SeqIO.write(seq_annot, salida, "genbank")


def parallel_runs(strains, loci_collection, threads):
    pool = multiprocessing.Pool(processes=threads)
    temp = partial(create_annotated_gbks, loci_collection=loci_collection)
    results = pool.map(func=temp, iterable=strains)
    pool.close()
    pool.join()


def tdna_classification(gbks_list, prefix, threads):
    loci_collection = identifica_locus(gbks_list)
    ant_to_analyze = set_anticodons(loci_collection)
    clustering_loci(ant_to_analyze, loci_collection)
    determine_nomenclature(loci_collection, ant_to_analyze, prefix)
    # create_annotated_gbks(gbks_list, loci_collection)
    parallel_runs(gbks_list, loci_collection, threads)