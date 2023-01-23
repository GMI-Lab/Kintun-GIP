#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 18:50:01 2020
Uses JolyTree to calculate a tree
@author: camilo
"""

import sys
import os
from distutils.spawn import find_executable
import logging


def is_tool(name):
    """Check whether `name` is on PATH."""
    return find_executable(name) is not None

def install_synerclust(path_file, path_synerclust):
    if os.path.isfile(path_synerclust):
        logging.info("Checking SynerClust: it was correctly detected")
    else:
        logging.info("Checking SynerClust: it was not detected")
        logging.info("Let's try to install SnerClust")
        install_path = os.path.join(os.path.abspath(path_file),
                                   'src/SynerClust/INSTALL.py')
        cmd_install = "python %s" % install_path
        os.popen(cmd_install).read()
        if os.path.isfile(path_synerclust):
            logging.info("Checking SynerClust: it was correctly detected")
        else:
            logging.info("SynerClust couldn't be installed. We're really sorry :(")
            sys.exit("SynerClust is not available for this code")
        

def run_jolytree(fasta_dir, dir_out, threads):
    # First check if jolytree is available
    for software in ['JolyTree.sh']:
        if is_tool(software):
            logging.info("""Checking %s: it was correctly detected
                         Running JolyTree with your fasta files"""
                         % software)
        else:
            logging.error("""Checking %s: ERROR: %s not found"""
                          % (software, software))
            sys.exit("""Checking %s: ERROR: %s not found"""
                     % (software, software))
    cmd1 = "cd %s; JolyTree.sh -i %s -b kintun-vli_tree -t %s; cd .." % (
            os.path.abspath(dir_out), os.path.abspath(fasta_dir), threads)
    os.popen(cmd1).read()


def run_synerclust(dir_out, threads, path_file):
    synerclust_path = os.path.join(os.path.abspath(path_file),
                                   'src/SynerClust/bin/synerclust.py')
    install_synerclust(path_file, synerclust_path)
    tree_path = os.path.join(os.path.abspath(dir_out), "kintun-vli_tree.nwk")
    syner_out = os.path.join(os.path.abspath(dir_out), "synerclust_output/")
    os.mkdir(syner_out)
    catalog_path = os.path.join(dir_out, "data_catalogue.txt")
    # -t tree.nwk [-n number_of_cores] [--run single]
    cmd2 = "%s -r %s -t %s -w %s -n %s --run single" % (
            synerclust_path, catalog_path, tree_path, syner_out, threads)
    os.popen(cmd2).read()
