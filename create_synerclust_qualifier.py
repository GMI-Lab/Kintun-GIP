#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 16:12:43 2020

@author: camilo
"""
from Bio import SeqIO


def create_dct(SC_file, number_strains):
    dct_clusters = {}
    core_genes = []
    with open(SC_file, 'r') as SC_results:
        for line in SC_results:
	    sep = line.strip().split("\t")
            if "(taxa: %s" % (str(number_strains)) in line:
                for i in sep[1].split():
                    dct_clusters[i.replace("_gene","")] = sep[0].split()[0]
                    core_genes.append(i.replace("_gene", ""))
    return dct_clusters, core_genes


def anotate_cluster_groups(gbks_list, SC_file):
    dctcluster, core_list = create_dct(SC_file, len(gbks_list))
    for strain_gbk in gbks_list:
        record = SeqIO.read(strain_gbk, 'genbank')
        new_features = []
        for feature in record.features:
            if feature.type == "CDS":
                if feature.qualifiers["locus_tag"][0] in dctcluster.keys():
                     key = feature.qualifiers["locus_tag"][0]
                     feature.qualifiers["synerclust_group"] = dctcluster[key]
                     new_features.append(feature)
                else:
                     continue
            elif feature.type == "gene":
                continue
            else:
                new_features.append(feature)
        record.features = new_features
        with open(strain_gbk.replace(".gbk","_SC_results.gbk"), "w") as salida:
            SeqIO.write(record, salida, "genbank")
