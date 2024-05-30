# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os
import pandas as pd

work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')


best_comb = 'P_1.0_C_1.0_G_0.0_M_10.0_E_10.0'


copy_from = 'step1_NMTF/output'
copy_to = f'{work_dir}/input'         
cell_conds = ['ND_8', 'ND_18', 'ND_25', 'ND_32', 'ND_37', 'WT_8', 'WT_18', 'WT_25', 'WT_32', 'WT_37']


for cc in cell_conds:
    for k1k2 in os.listdir(f'{copy_from}/{cc}'):
        for file in os.listdir(f'{copy_from}/{cc}/{k1k2}/{best_comb}'):
            if file == 'G1_with_headers.csv' or file == 'S_EXP.csv':
                # case = root.split('/')[5]
                sd = f'{copy_to}/NMTF_G1s/{best_comb}/{cc}'
                if not os.path.exists(sd):
                    os.makedirs(sd)
                copyfile(f'{copy_from}/{cc}/{k1k2}/{best_comb}/{file}', f'{sd}/{cc}_{file}')   

# for root, dirs, files in os.walk(copy_from):
#     for file in files:
#         if 'G1_with_headers' in file or 'S_EXP' in file:
#             if 'ALL' in root:
#                 cell_cond = root.split('\\')[1]
#                 if not os.path.exists( f'{copy_to}/{cell_cond}'):
#                     os.makedirs(f'{copy_to}/{cell_cond}')
#                 copyfile(f'{root}/{file}', f'{copy_to}/{cell_cond}/{file}')


os.chdir(work_dir) 
