# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 16:57:53 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os


work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

copy_from = 'step2_StageSpecPreds/output'
copy_to = f'{work_dir}/input'  

for file in os.listdir(copy_from):
    if '.pkl' in file and 'Genes' in file:
        if not os.path.exists(f'{copy_to}/Geneslist'):
            os.makedirs(f'{copy_to}/Geneslist')
        copyfile(f'{copy_from}/{file}', f'{copy_to}/Geneslist/{file}')
    elif '.pkl' in file and 'Genes' not in file:
        copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')
       


path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


copyfile('Data/LitValid_AllGenes/output/LitValid_AllGenes.pkl', f'{work_dir}/input/LitValid_AllGenes.pkl') 
# copyfile('Data/PDmap/output/PDgenes_PDmap_noUK.pkl', f'{work_dir}/input/PDgenes_PDmap_noUK.pkl') 
# copyfile('Data/PDmap/output/PDmap_BasePathways.pkl', f'{work_dir}/input/PDmap_BasePathways.pkl') 
# copyfile('Data/PDmap/output/PDgenes_PDmap.pkl', f'{work_dir}/input/PDgenes_PDmap.pkl') 
copyfile('Data/Parse_Gene4PD/output/Gene4PD.pkl', f'{work_dir}/input/Gene4PD.pkl') 

copyfile('Data/Homo_sapiens.gene_info', f'{work_dir}/input/Homo_sapiens.gene_info') 

os.chdir(work_dir)



