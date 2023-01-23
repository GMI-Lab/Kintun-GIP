#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:18:22 2020
Write catalogue for Synerclust
@author: camilo
"""

import os
import ntpath
from Bio import SeqIO
import logging


def create_catalogue(file_out):
    if os.path.isfile(file_out):
        logging.info("data_catalogue.txt was overwritten")
        logging.info("---------------")
    with open(file_out, 'w') as out_file:
        out_file.write("transl_table\t11\n")
        out_file.write("//\n")


def genomes_to_catalogue(list_fasta, dir_out, ext):
    out_file = os.path.join(dir_out, "data_catalogue.txt")
    create_catalogue(out_file)
    for file in list_fasta:
        prefix = ntpath.basename(file).split(".%s" % ext)[0]
        with open(out_file, 'a') as cata_file:
            # Write strain code based in filename in catalogue
            cata_file.write("Genome\t%s\n" % prefix)
            # Write the genomic fasta file
            prokka_dir = os.path.join(
                dir_out, ntpath.basename(file)+"_prokka/")
            cata_file.write("Sequence\t%s\n" % os.path.abspath(
                            os.path.join(prokka_dir, prefix + ".fna")))
            # Write the annotation gff file
            cata_file.write("Annotation\t%s\n" % os.path.abspath(
                            os.path.join(prokka_dir, prefix + ".gff")))
            cata_file.write("//\n")
    # Write the last line of the catalogue
    with open(out_file, 'a') as cata_file:
        cata_file.write("//\n")
    logging.info("data_catalogue.txt was created")
    logging.info("---------------")
