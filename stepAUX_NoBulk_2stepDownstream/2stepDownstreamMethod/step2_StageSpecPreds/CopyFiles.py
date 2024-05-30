# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 11:34:01 2023

@author: kmihajlo
"""


from shutil import copyfile, copytree, rmtree
import os

def split_file_path(file_path):
    # Normalize the file path to use the appropriate separator for the current OS
    normalized_path = os.path.normpath(file_path)
    # Split the file path into its components
    components = normalized_path.split(os.sep)
    return components


work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')
 
    
copy_from = 'step1_PDEnrClstsInitPreds/output'  
copy_to = f'{work_dir}/input/InitPreds'  
if not os.path.exists(copy_to):
    os.makedirs(copy_to)

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if 'InitPreds.pkl' in file and 'ND' in file:
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')
            
    
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
 
copy_from = 'step2_NMTF_WeightedNets/output'
copy_to = f'{work_dir}/input/NMTF_matrices'         
best_comb = 'P_1.0_C_1.0_G_0.0_M_10.0_E_10.0'     

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file == 'G1_with_headers.csv' and best_comb in root:
            cc = split_file_path(root)[2]
            if not os.path.exists(f'{copy_to}/{cc}'):
                os.makedirs(f'{copy_to}/{cc}')

            copyfile(f'{root}/{file}', f'{copy_to}/{cc}/{cc}_{file}')
        if file == 'S_EXP.csv' and best_comb in root:
            copyfile(f'{root}/{file}', f'{copy_to}/{cc}/{cc}_{file}')

copyfile('Data/DisGeNet/output/PDgenes_DGN_ALL.pkl', f'{work_dir}/input/PDgenes_DGN_ALL.pkl') 

copy_from = 'step7_TGM_Predictions/output/D8_D18_D25_D32_D37/1.533perc'
copy_to = f'{work_dir}/input'  
file = 'GM_D8_D18_D25_D32_D37_1.533_PDpreds.csv'
copyfile(f'{copy_from}/{file}', f'{copy_to}/MO_{file}')



os.chdir(work_dir) 
