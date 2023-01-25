#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 18:34:12 2020
This script run PROKKA over a genome, and modify the gff file
for Synerclust
@author: camilo
"""
import sys
import os
from Bio import SeqIO
from distutils.spawn import find_executable
import logging
import ntpath
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
#from Bio.Alphabet import IUPAC

import shutil


def is_tool(name):
    """Check whether `name` is on PATH."""
    return find_executable(name) is not None


def annotate_genome(list_fasta, dir_out, file_ext, threads):
    # First check if Prodigal, aragorn and barrnap is available
    for software in ['aragorn', 'barrnap']:
        if is_tool(software):
            logging.info("""Checking %s: it was correctly detected"""
                         % software)
        else:
            logging.error("""Checking %s: ERROR: %s not found"""
                          % (software, software))
            sys.exit("""Checking %s: ERROR: %s not found"""
                     % (software, software))
    # Run softwares
    for file in list_fasta:
        prefix = ntpath.basename(file).split(".%s" % file_ext)[0]
        outdir = os.path.join(dir_out, ntpath.basename(file)+"_annot/")
        os.mkdir(outdir)
        logging.info("Processing %s file %s/%s" % (ntpath.basename(file),
                     str(list_fasta.index(file)+1), str(len(list_fasta))))

        cmd2 = "aragorn -fo -m -t -gcbact -c -fasta %s > %s/%s.prev" %(file, outdir, prefix)
        os.popen(cmd2).read()
        # Process stdout from Aragorn
        counter = 1
        out_file = "%s/%s.aragorn" %(outdir, prefix)
        with open("%s/%s.prev" %(outdir, prefix), "r") as trnas_in:
            for line in trnas_in:
                if line[0] == ">":
                    sep = line.strip().split(" ")
                    tdna_id = sep[0].replace(">", "")
                    tdna_start = sep[1].split("[")[1].split(",")[0]
                    tdna_end = sep[1].split(",")[1].split("]")[0]
                    if sep[1][0] == "c":
                        tdna_sense = "-"
                    else:
                        tdna_sense = "+"
                    line1 = "%s\tARAGORN\tt(m)DNA\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                    prefix, tdna_start, tdna_end, tdna_sense, str(counter).zfill(4), tdna_id)
                    counter += 1
                    with open(out_file, 'a') as salida:
                        salida.write(line1)
                else:
                    continue
        os.remove("%s/%s.prev" %(outdir, prefix))
        # Run barrnap
        cmd3 = "barrnap --kingdom bact --threads %s --quiet %s > %s/%s.barrnap" %(threads, file, outdir, prefix)
        os.popen(cmd3).read()

