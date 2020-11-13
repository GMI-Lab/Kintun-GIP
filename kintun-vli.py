#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 14:59:58 2020

@author: camilo
"""

import os
import sys
import argparse
import logging
import time
import check_fasta
import annotate_fasta
import shutil
import create_data_catalogue
import run_synerclust
import glob2
import create_synerclust_qualifier
import tdna_classification

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, action="store", dest='fasta_dir',
                        help="Directory with fasta files (REQUIRED)", required=True)
    parser.add_argument("-x", type=str, action="store", dest='file_extension',
                        help="Extension for fasta files (REQUIRED;default: fa)",
                        default='fa', required=True)
    parser.add_argument("-o", type=str, action="store", dest='results_dir',
                        default='', required=True,
                        help="Name for results directory (REQUIRED)")
    parser.add_argument("-p", type=str, action="store", dest='prefix',
                        default='NNN', required=True,
                        help="Prefix for t(m)DNA scheme classification (REQUIRED)")
    parser.add_argument("-t", type=str, action="store", dest='threads',
                        default='2',
                        help="Number of threads for PROKKA and SYNERCLUST")
    parser.add_argument("-log", type=str, action="store", dest='log_file',
                        default='',
                        help="Name of log file (Optional)")
    parser.add_argument("-v", action="version",
                        version='Kintun-VLI tDNAs-classification v1.0')
    args = parser.parse_args()
    if os.path.exists(args.results_dir):
        sys.exit("Output directory already exists!")
    os.makedirs(args.results_dir)
    if args.log_file == '':
        log_file = time.strftime("kintun-vli_%Y%m%d-%H%M.log")
    else:
        log_file = args.log_file
    logging.basicConfig(filename=os.path.join(args.results_dir, log_file),
                        level=20)
    logging.info('This is Kintun-VLI for tDNAs and tmDNAs classification')
    logging.info('If you use it, please cite us')
    logging.info("---------------")
    logging.info("You're awesome!! :)")
    logging.info("---------------")
    logging.info("Let's create the output directory")
    logging.info("%s was created" % args.results_dir)
    logging.info("Let's check your fasta files")
    logging.info("---------------")
    # Check fasta files
    list_files = [os.path.join(args.fasta_dir, file)
                  for file in os.listdir(args.fasta_dir)
                  if file.endswith(args.file_extension)]
    if len(list_files) < 4:
        logging.info("Not enough fasta files to analyze detected! Exit!")
        sys.exit("Not enough fasta files to analyze detected! Exit!")
    else:
        logging.info("Were detected %s fasta files" % str(len(list_files)))
    for file in list_files:
        if check_fasta.is_fasta(file):
            logging.info("%s looks good!" % file)
        else:
            logging.error(
                "Please check file with name %s is a fasta file!" % file)
            sys.exit("There is something wrong with your fasta files :(")
    logging.info("---------------")
    logging.info("Now let's annotate your genomes with PROKKA")
    logging.info("---------------")

    annotate_fasta.annotate_genome(list_files, args.results_dir,
                                   args.file_extension, args.threads)
    logging.info("---------------")
    logging.info("All your genomes are now annotated, yay! :D")
    logging.info("Let's fix some things for Synerclust to run")
    logging.info("---------------")
    create_data_catalogue.genomes_to_catalogue(list_files, args.results_dir,
                                               args.file_extension)
    run_synerclust.run_jolytree(args.fasta_dir, args.results_dir, args.threads)
    logging.info("---------------")
    logging.info("All seems good, let's run synerclust!")
    logging.info("---------------")
    run_synerclust.run_synerclust(args.results_dir, args.threads,
                                  os.path.abspath(__file__).replace(
                                      "kintun-vli.py", ""))
    logging.info("---------------")
    logging.info("Finally it's time for tDNAs classification!")
    logging.info("First, Synerclust groups have to be added to gbks")
    logging.info("---------------")
    create_synerclust_qualifier.anotate_cluster_groups(
        glob2.glob(os.path.abspath(args.results_dir) + "/*prokka/*.gbk", recursive=True), os.path.join(
                  os.path.abspath(args.results_dir),
                  "synerclust_output/results/final_clusters.txt"))
    logging.info("---------------")
    logging.info("And now, tDNAs classification!")
    logging.info("If your FASTA collection is big it can take a while")
    logging.info("---------------")
    tdna_classification.tdna_classification(
        glob2.glob(os.path.abspath(args.results_dir) +
                  "/*prokka/*_SC_results.gbk", recursive=True), args.prefix)
    logging.info("---------------")
    logging.info("It seems classification is over, you can find tDNAS annotated in the final genbank files")
    logging.info("Have a good day!")
    logging.info("---------------")

if __name__ == '__main__':
    main()
