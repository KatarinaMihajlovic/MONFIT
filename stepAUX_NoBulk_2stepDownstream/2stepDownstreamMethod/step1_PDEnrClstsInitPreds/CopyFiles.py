# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:49:13 2023

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os

work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

copy_from = 'step3_ClusterG1NMTF_WeightedNets/output'
best_comb = 'P_1.0_C_1.0_G_0.0_M_10.0_E_10.0'     

for cc in os.listdir(copy_from):
    copy_to = f'{work_dir}/input/Clusters/{cc}' 
    if not os.path.exists(copy_to):
        os.makedirs(copy_to)   
    for file in os.listdir(f'{copy_from}/{cc}/{best_comb}'):
        copyfile(f'{copy_from}/{cc}/{best_comb}/{file}', f'{copy_to}/{file}')

    
copyfile('Data/DisGeNet/output/PDgenes_DGN_ALL.pkl', f'{work_dir}/input/PDgenes_DGN_ALL.pkl') 
os.chdir(work_dir) 
