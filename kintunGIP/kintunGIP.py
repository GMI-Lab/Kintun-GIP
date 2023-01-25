#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.append('/')

import os
import argparse
import logging
import time
import check_fasta
import annotate_fasta
import tdna_clusterization

class Kintun:
    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Kintun-GIP is a genome islands predictor on bacterial genomes',
            usage='''kintunGIP.py <command> [<args>]

    The available commands are:
      clust     Clusterize and name all t(m)DNAs in a bacterial genome collection.
       (Coming soon)
    ''')

        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    @staticmethod
    def clust():
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", type=str, action="store", dest='fasta_dir',
                            help="Directory with fasta files (REQUIRED)", required=True)
        parser.add_argument("-x", type=str, action="store", dest='file_extension',
                            help="Extension for fasta files (REQUIRED;default: fasta)",
                            default='fasta', required=True)
        parser.add_argument("-o", type=str, action="store", dest='results_dir',
                            default='', required=True,
                            help="Name for results directory (REQUIRED)")
        parser.add_argument("-p", type=str, action="store", dest='prefix',
                            default='NNN', required=True,
                            help="Prefix for the t(m)DNA clustering scheme (REQUIRED)")
        parser.add_argument("-t", type=str, action="store", dest='threads',
                            default='8',
                            help="Number of threads for nucmer and SibeliaZ")
        parser.add_argument("-log", type=str, action="store", dest='log_file',
                            default='',
                            help="Name of log file (Optional)")
        parser.add_argument("-v", action="version",
                            version='Kintun-GIP t(m)DNAs clusterization module v1.0')
        args = parser.parse_args(sys.argv[2:])
        if os.path.exists(args.results_dir):
            sys.exit("Output directory already exists!")
        os.makedirs(args.results_dir)

        if args.log_file == '':
            log_file = time.strftime("kintun-vli_%Y%m%d-%H%M.log")
        else:
            log_file = args.log_file

        logging.basicConfig(filename=os.path.join(args.results_dir, log_file),
                            level=20)
        logging.info('This is Kintun-VLI for tDNAs and tmDNAs clusterization')
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
        for file_ii in list_files:
            if check_fasta.is_fasta(file_ii):
                logging.info("%s looks good!" % file_ii)
            else:
                logging.error(
                    "Please check if file %s is a fasta file!" % file_ii)
                sys.exit("There is something wrong with one of your fasta files :(")
        logging.info("---------------")
        logging.info("Now let's annotate your genomes with Aragorn and Barrnap")
        logging.info("---------------")
        annotate_fasta.annotate_genome(list_files, args.results_dir,
                                       args.file_extension, args.threads)
        logging.info("---------------")
        logging.info("All your genomes are now annotated, yay! :D")
        logging.info("Finally it's time for tDNAs clusterization")
        tdna_clusterization.tdna_clusterization(args.fasta_dir, args.results_dir, args.file_extension, args.prefix)
        #logging.info("---------------")
        #logging.info("It seems clusterization is over, you can find tDNAS annotated in the final genbank files")
        #logging.info("Have a good day!")
        #logging.info("---------------")


if __name__ == '__main__':
    Kintun()
