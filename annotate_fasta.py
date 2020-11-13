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
import shutil


def is_tool(name):
    """Check whether `name` is on PATH."""
    return find_executable(name) is not None


def annotate_genome(list_fasta, dir_out, file_ext, threads):
    # First check if PROKKA is available
    for software in ['prokka']:
        if is_tool(software):
            logging.info("""Checking %s: it was correctly detected"""
                         % software)
        else:
            logging.error("""Checking %s: ERROR: %s not found"""
                          % (software, software))
            sys.exit("""Checking %s: ERROR: %s not found"""
                     % (software, software))
    # Run PROKKA
    for file in list_fasta:
        prefix = ntpath.basename(file).split(".%s" % file_ext)[0]
        outdir = os.path.join(dir_out, ntpath.basename(file)+"_prokka/")
        logging.info("Processing %s file %s/%s" % (ntpath.basename(file),
                     str(list_fasta.index(file)+1), str(len(list_fasta))))
        cmd1 = "prokka --addgenes --cpus %s --fast --quiet --outdir %s --prefix %s %s" % (
            threads, outdir, prefix, file)
        os.popen(cmd1).read()
        logging.info("%s processed!" % ntpath.basename(file))
        for result in os.listdir(outdir):
            # Correct gff file for Synerclust
            if os.path.splitext(result)[1] == ".gff":
                with open(outdir + result, 'r') as in_file:
                    with open(outdir + result + ".modified", 'w') as out_file:
                        for line in in_file:
                            sep = line.strip().split("\t")
                            if len(sep) > 1:
                                if sep[2] == 'gene':
                                    out_file.write(line)
                                elif sep[2] == 'CDS':
                                    out_file.write(line)
                                    out_file.write(line.replace("CDS", "mRNA"))
                                    out_file.write(line.replace("CDS", "exon"))
                                    out_file.write("\n")
                                else:
                                    continue
                shutil.move(outdir + result + ".modified", outdir + result)
            elif os.path.splitext(result)[1] in [".gbk", ".fna"]:
                continue
            # remove unused files
            else:
                os.remove(os.path.join(outdir, result))
