# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:24:44 2021

@author: roberto
"""

import sys
import subprocess
import csv
import pandas as pd
from Bio import SeqIO
from Bio import SeqFeature as sf
from Bio.SeqRecord import SeqRecord
import glob
import os
import datetime
import json

# Primera función, chequear las coordenadas

def extract_pilrs(original_pilrs_df):
    # Create a new DataFrame
    new_data = {
        'strain': [],
        'start': [],
        'end': [],
        'sense': [],
        'tdna_name': [],
        'condition_type': [],
        'condition_start': [],
        'condition_end': [],
    }
    # Open the PILRs dataframe
    df = pd.read_csv(original_pilrs_df, header=True, index=0)

    # Iterate rows in df
    for index, row in df.iterrows():
        if row['ur_condition'] != 'Empty':
            new_data['strain'].append(row['strain'])
            new_data['start'].append(row['start'])
            new_data['end'].append(row['end'])
            new_data['sense'].append(row['sense'])
            new_data['tdna_name'].append(row['tdna_name'])
            new_data['condition_type'].append('ur_condition')
        if row['dr_condition'] != 'Empty':
            new_data['strain'].append(row['strain'])
            new_data['start'].append(row['start'])
            new_data['end'].append(row['end'])
            new_data['sense'].append(row['sense'])
            new_data['tdna_name'].append(row['tdna_name'])
            new_data['condition_type'].append('dr_condition')

    # Return the new DataFrame
    return pd.DataFrame(new_data)

def pilrs_annotation(input_folder, threads, prefix):
    # Step 1: Run SibeliaZ to identify LCBs and create synteny blocks and create LCBS dataframe
    file_pilrs = glob.glob(f'{input_folder}/*tmDNAS_PILRs.csv')
    pilrs_df = extract_pilrs(file_pilrs[0])