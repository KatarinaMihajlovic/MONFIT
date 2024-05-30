# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:25:21 2024

@author: Katarina
"""

from shutil import copyfile
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


copy_from = 'step1_CreateNetworks/output'

copy_to = f'{work_dir}/input'

cell_conds = ['ND_18', 'WT_18', 'WT_32', 'ND_32', 'WT_8', 'ND_8', 'ND_25', 'WT_25', 'WT_37', 'ND_37']
for cc in cell_conds:
    save_dir = f'{copy_to}/{cc}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    copyfile(f'{copy_from}/{cc}/E_{cc}_normalized.csv', f'{save_dir}/E_{cc}_normalized.csv')

copyfile('Data/Homo_sapiens.gene_info', f'{work_dir}/input/Homo_sapiens.gene_info')

copy_from = 'Data/MolecularNetworks'
copy_to = f'{work_dir}/input/MolecularNetworks'
if not os.path.exists(copy_to):
    os.makedirs(copy_to)
    
copyfile(f'{copy_from}/output/Human_PPI_General.edgelist', f'{copy_to}/Human_PPI_General.edgelist')
copyfile(f'{copy_from}/output/Human_COEX_General.edgelist', f'{copy_to}/Human_COEX_General.edgelist')
copyfile(f'{copy_from}/output/Human_GI_General.edgelist', f'{copy_to}/Human_GI_General.edgelist')
copyfile(f'{copy_from}/output/enzyme_subprodgene_KEGG.csv', f'{copy_to}/enzyme_subprodgene_KEGG.csv')


os.chdir(work_dir) 
