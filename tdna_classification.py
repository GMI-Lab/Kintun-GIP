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


class locus():
    def __init__(locus, strain, tag):
        locus.strain = strain
        locus.tag = "%s_%s" % (strain, tag)
        locus.sense = ''
        locus.tdnas = []
        locus.anticodon_in_locus = []
        locus.UR_posI = ''
        locus.UR_posII = ''
        locus.UR_posIII = ''
        locus.UR_groups = []
        locus.pos_senses = []
        locus.tdnas_in_UR = ''
        locus.rdnas_in_UR = ''
        locus.UR_I_class = []
        locus.UR_II_class = []
        locus.UR_III_class = []
        locus.classification = []
        locus.context_list = []
        locus.record = []


def identifica_locus(strain_codes):
    dct_ant, aa_cod = create_dics()
    loci_collection = []
    tag = 1
    for strain in strain_codes:
        seq_annot = SeqIO.read(strain, 'genbank')
        # Crea las listas de anotaciones
        annotations_list = seq_annot.features
        # Busca en las listas los tDNAs o tmDNAs (sí, también tmDNAs, uhlalá)
        x = 0
        while x < len(annotations_list) - 1:
            feature = annotations_list[x]
            if feature.type in ['tRNA', 'tmRNA']:
                new_locus = locus(strain, tag)
                new_locus.sense = feature.strand
                new_locus.tdnas.append(feature)
                new_locus.record = feature
                for idx in range(x+1, len(annotations_list)):
                    if annotations_list[idx].type == 'tRNA' and new_locus.sense == annotations_list[idx].strand:
                        new_locus.tdnas.append(annotations_list[idx])
                        continue
                    else:
                        break
                # Coordenadas del locus
                locus_end = idx
                x = idx
                # Identifica si es un operon o una unidad mono o policistrónica
                if len([new_locus.tdnas]) > 1:
                    new_locus.classification = 'poly-tu'
                else:
                    new_locus.classification = 'mono-tu'
                # Identifica las tres posiciones río arriba
                UR_sequences = []
                tdnas_UR = []
                rdnas_UR = []
                for feature in (annotations_list[:x][::-1] +
                                annotations_list[::-1]):
                    if len(UR_sequences) < 3:
                        if feature.type == 'CDS':
                            UR_sequences.append(feature)
                        elif feature.type == 'tRNA':
                            tdnas_UR.append(
                                aa_cod[feature.qualifiers['product'][0]])
                        elif feature.type == 'rRNA':
                            rdnas_UR.append(feature.qualifiers['product'][0])
                        else:
                            continue
                    else:
                        break
                DR_sequences = []
                tdnas_DR = []
                rdnas_DR = []
                for feature in annotations_list[idx:] + annotations_list:
                    if len(DR_sequences) < 3:
                        if feature.type == 'CDS':
                            DR_sequences.append(feature)
                        elif feature.type == 'tRNA':
                            tdnas_DR.append(
                                aa_cod[feature.qualifiers['product'][0]])
                        elif feature.type == 'rRNA':
                            rdnas_DR.append(feature.qualifiers['product'][0])
                        else:
                            continue
                    else:
                        break
                posi_keys = ['I', 'II', 'III']
                for k in [0, 1, 2]:
                    if new_locus.sense == 1:
                        exec(
                            "new_locus.UR_pos%s = UR_sequences[%s]" %
                            (posi_keys[k], str(k))
                            )
                        new_locus.tdnas_in_UR = tdnas_UR
                        new_locus.rdnas_in_UR = rdnas_UR
                    else:
                        exec(
                            "new_locus.UR_pos%s = DR_sequences[%s]" %
                            (posi_keys[k], str(k))
                            )
                        new_locus.tdnas_in_UR = tdnas_DR
                        new_locus.rdnas_in_UR = rdnas_DR
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
                x = locus_end
                tag = tag + 1
            else:
                x = x + 1
                continue
    return loci_collection


def set_anticodons(loci_collection):
    ants_in_loci = []
    for tdna in loci_collection:
        for i in tdna.anticodon_in_locus:
            ants_in_loci.append(i)
    ant_out = list(set(ants_in_loci))
    return ant_out


def generate_binary_data(anticodon, loci_to_analyze):
    posI_class = []
    posII_class = []
    posIII_class = []
    for loci in loci_to_analyze:
        posI_class.append(loci.UR_groups[0])
        posII_class.append(loci.UR_groups[1])
        posIII_class.append(loci.UR_groups[2])
    list_categories = list(set(posI_class + posII_class + posIII_class))
    binary_pos_collection = []
    for loci in loci_to_analyze:
        binary_loci = []
        classes_in_loci = loci.UR_groups
        for cat in list_categories:
            if cat in classes_in_loci:
                binary_loci.append(1)
            else:
                binary_loci.append(0)
        binary_pos_collection.append(binary_loci)
    return binary_pos_collection


def clustering_loci(ant_to_analyze, all_loci):
    dct_ant, aa_cod = create_dics()
    # Select loci
    for ant_cod in ant_to_analyze:
        coll_loci = []
        for loci in all_loci:
            if ant_cod in loci.anticodon_in_locus:
                coll_loci.append(loci)
            else:
                continue
        if len(coll_loci) > 1:
            sequences_by_strain = {}
            for loci in coll_loci:
                sequences_by_strain[loci.strain] = []
            binary_data = generate_binary_data(ant_cod, coll_loci)
            labels = fclusterdata(binary_data, len(coll_loci),
                                  criterion='maxclust', metric='hamming',
                                  method='ward')
            for index in range(len(coll_loci)):
                coll_loci[index].context_list.append("%s_%s" % (
                    dct_ant[ant_cod], labels[index]))
        else:
            coll_loci[0].context_list.append("%s_%s" % (dct_ant[ant_cod], '1'))


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
        coll_loci = []
        for locus in all_loci:
            if ant_cod in locus.anticodon_in_locus:
                coll_loci.append(locus)
        contexts_types = []
        for locus in coll_loci:
            for i in locus.context_list:
                if dct_ant[ant_cod] in i and i not in contexts_types:
                    contexts_types.append(i)
                else:
                    continue
        contexts_types_info = []
        for i in contexts_types:
            strains = []
            list_number_tdnas = []
            for locus in coll_loci:
                if i in locus.context_list:
                    strains.append(locus.strain)
                    number_tdnas = []
                    for element in locus.tdnas:
                        if element.type == "tRNA" and element.qualifiers['note'][0].split('(')[1].strip(')').upper() == ant_cod:
                            number_tdnas.append(element)
                        elif element.type == "tmRNA":
                            number_tdnas.append(element)
                        else:
                            continue
                list_number_tdnas.append(len(number_tdnas))
            contexts_type = context_type(i, len(strains), np.max(list_number_tdnas))
            contexts_types_info.append(contexts_type)
        # 1
        # 2
        # 3
        contexts_types_info = sorted(contexts_types_info, key=lambda context_type:context_type.max_tdnas, reverse=True)
        contexts_types_info = sorted(contexts_types_info, key=lambda context_type:context_type.num_strains, reverse=True)
        list_alpha = create_letter_list()
        counter = 0
        for context in contexts_types_info:
            for i in range(counter, counter + context.max_tdnas):
                context.letters.append(list_alpha[i])
            counter = counter + context.max_tdnas
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


def create_annotated_gbks(strain_codes, loci_collection):
    for strain in strain_codes:
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


def tdna_classification(gbks_list, prefix):
    loci_collection = identifica_locus(gbks_list)
    ant_to_analyze = set_anticodons(loci_collection)
    clustering_loci(ant_to_analyze, loci_collection)
    determine_nomenclature(loci_collection, ant_to_analyze, prefix)
    create_annotated_gbks(gbks_list, loci_collection)
