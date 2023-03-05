#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 18:34:12 2020
This script run tRNAscan-SE and BARRNAP over a chromosome
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
from tqdm.auto import tqdm
#from Bio.Alphabet import IUPAC

import shutil


def is_tool(name):
    """Check whether `name` is on PATH."""
    return find_executable(name) is not None


def annotate_genome(list_fasta, dir_out, file_ext, threads):
    # First check if Prodigal, aragorn and barrnap is available
    for software in ['tRNAscan-SE', 'aragorn', 'barrnap']:
        if is_tool(software):
            logging.info("""Checking %s: it was correctly detected"""
                         % software)
        else:
            logging.error("""Checking %s: ERROR: %s not found"""
                          % (software, software))
            sys.exit("""Checking %s: ERROR: %s not found"""
                     % (software, software))
    # Run softwares
    for in_file in tqdm(list_fasta):
        prefix = ntpath.basename(in_file).replace(f".{file_ext}", "")
        outdir = os.path.join(dir_out, ntpath.basename(in_file).replace(f".{file_ext}", "")+"_annot/")
        new_fasta = f"{outdir}/{prefix}.fasta"
        os.mkdir(outdir)
        #Create new fasta, avoids input format problems with tRNAscan-SE
        with open(in_file,"r") as ori_fasta:
            for record in SeqIO.parse(ori_fasta,"fasta"):
                with open(new_fasta,"w") as fasta_process:
                    SeqIO.write(record,fasta_process,"fasta")
        logging.info("Processing %s file %s/%s" % (ntpath.basename(in_file),
                     str(list_fasta.index(in_file)+1), str(len(list_fasta))))
        # Run tRNAscan-SE
        cmd2 = f"tRNAscan-SE -qQ -B {new_fasta} > {outdir}/{prefix}.prev"
        os.popen(cmd2).read()
        # Process stdout from tRNAscan-SE
        counter = 1
        out_file = "%s/%s.tmdnas" %(outdir, prefix)
        with open("%s/%s.prev" %(outdir, prefix), "r") as trnas_in:
            for line in trnas_in:
                if "Sequence" in line or "Name" in line or "--------" in line or "NNN" in line:
                    continue
                else:
                    sep = line.strip().split("\t")
                    if int(sep[2]) < int(sep[3]):
                        tdna_sense = "+"
                        tdna_start = sep[2]
                        tdna_end = sep[3]
                    else:
                        tdna_sense = "-"
                        tdna_start = sep[3]
                        tdna_end = sep[2]
                    #(0)Kp_1084 	1	221662 	221738 	Pro	CGG	0	0	76.0
                    line1 = "%s\ttRNAscan-SE\tt(m)DNA\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s-%s\n" % (
                    prefix, tdna_start, tdna_end, tdna_sense, sep[1].zfill(4), sep[4], sep[5])
                    counter += 1
                    with open(out_file, 'a') as salida:
                        salida.write(line1)
        os.remove("%s/%s.prev" %(outdir, prefix))
        #Run ARAGORN
        cmd3 = f"aragorn -fo -m -gcbact -c -fasta {new_fasta} > {outdir}/{prefix}.prev"
        os.popen(cmd3).read()
        # Process stdout from Aragorn
        counter = 1
        with open("%s/%s.prev" % (outdir, prefix), "r") as trnas_in:
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
        os.remove("%s/%s.prev" % (outdir, prefix))
        # Run barrnap
        cmd3 = "barrnap --kingdom bact --threads %s --quiet %s > %s/%s.rrnas" %(threads, in_file, outdir, prefix)
        os.popen(cmd3).read()


