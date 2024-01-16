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
import ast

def extract_pilrs(original_pilrs_df):
    # Create a new DataFrame
    new_data = {
        'strain': [],
        'start': [],
        'end': [],
        'sense': [],
        'tdna_name': [],
        'pilr_type': [],
        'condition_type': [],
    }
    # Open the PILRs dataframe
    df = pd.read_csv(original_pilrs_df, header=0, index_col=0)

    # Iterate rows in df
    for index, row in df.iterrows():
        if row['ur_condition'] != 'Empty':
            new_data['strain'].append(row['strain'])
            new_data['start'].append(row['start'])
            new_data['end'].append(row['end'])
            new_data['sense'].append(row['sense'])
            new_data['tdna_name'].append(row['tdna_name'])
            new_data['pilr_type'].append('UR_PILR')
            new_data['condition_type'].append(row['ur_condition'])

        if row['dr_condition'] != 'Empty':
            new_data['strain'].append(row['strain'])
            new_data['start'].append(row['start'])
            new_data['end'].append(row['end'])
            new_data['sense'].append(row['sense'])
            new_data['tdna_name'].append(row['tdna_name'])
            new_data['pilr_type'].append('DR_PILR')
            new_data['condition_type'].append(row['dr_condition'])

    # Return the new DataFrame
    return pd.DataFrame(new_data)

def convert_to_tuple(string):
    return ast.literal_eval(string)


def process_overlapping_rows(group):
    result_rows = []

    for _, row in group.iterrows():
        overlapping_rows = []
        for _, other_row in group.iterrows():
            if check_if_overlap((row['PILR_start'],row['PILR_end']), (other_row['PILR_start'],other_row['PILR_end'])):
                overlapping_rows.append(other_row)
        if overlapping_rows:
            # Retain the row with the nearest PILR
            nearest_row = min(overlapping_rows, key=lambda r: abs(row['start'] - r['PILR_start']) if row['sense'] == '+' else abs(row['end'] - r['PILR_end']))
            result_rows.append(nearest_row)
        else:
            result_rows.append(row)

    return pd.DataFrame(result_rows)


def pilrs_annotation(input_folder, threads, prefix):
    # Step 1: Run SibeliaZ to identify LCBs and create synteny blocks and create LCBS dataframe
    file_pilrs = glob.glob(f'{input_folder}/*tmDNAS_PILRs.csv')
    pilrs_df = extract_pilrs(file_pilrs[0])
    pilrs_df[['PILR_start', 'PILR_end', 'PILR_length']] = pilrs_df['condition_type'].apply(convert_to_tuple).apply(
        pd.Series)
    pilrs_df = pilrs_df.drop('condition_type', axis=1)

    # Apply the function to each strain group
    pilrs_df = pilrs_df.groupby('strain').apply(process_overlapping_rows).reset_index(drop=True)

    pilrs_df = pilrs_df.drop_duplicates()

    pilrs_df.to_csv(f"{input_folder}/{prefix}_PILRs_details.csv")
