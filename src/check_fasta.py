#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:54:10 2020
Check if a given file is in FastA format
@author: CBP
"""
from Bio import SeqIO


def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty
