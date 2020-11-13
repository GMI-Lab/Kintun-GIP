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


def run_jolytree(fasta_dir, dir_out, threads):
    # First check if PROKKA is available
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
