# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""

from shutil import copyfile 
import os, sys

def split_file_path(file_path):
    # Normalize the file path to use the appropriate separator for the current OS
    normalized_path = os.path.normpath(file_path)
    # Split the file path into its components
    components = normalized_path.split(os.sep)
    return components

E_w = float(sys.argv[1])
PPI_w = float(sys.argv[2])
COEX_w = float(sys.argv[3])
GI_w = float(sys.argv[4])
MI_w = float(sys.argv[5])
wd = str(sys.argv[6])

if not os.path.exists(f'{wd}/input'):
    os.makedirs(f'{wd}/input') 
if not os.path.exists(f'{wd}/output'):
    os.makedirs(f'{wd}/output')

best_comb = f'P_{PPI_w}_C_{COEX_w}_G_{GI_w}_M_{MI_w}_E_{E_w}'


copy_from = 'step2_NMTF_WeightedNets/output'
copy_to = f'{wd}/input'         

for root, dirs, files in os.walk(copy_from):
    if best_comb in root:
        for file in files:
            if file == 'G1_with_headers.csv':
                cc = split_file_path(root)[2]

                # case = root.split('/')[5]
                sd = f'{copy_to}/NMTF_G1s/{cc}/{best_comb}'
                if not os.path.exists(sd):
                    os.makedirs(sd)
                copyfile(f'{root}/{file}', f'{sd}/{cc}_{file}')   


