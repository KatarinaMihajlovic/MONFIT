# -*- coding: utf-8 -*-
"""
Created on Thu Oct  14 12:48:56 2021

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


# copy_to = f'{work_dir}/input'                 
# dir2copy = 'step7_SimMeasuresG1S5/output'
# if os.path.exists(f'{copy_to}/SimMeasures'): 
#     rmtree(f'{copy_to}/SimMeasures') 
# copytree(dir2copy, f'{copy_to}/SimMeasures') 


# copyfile('CommonData/Gene4PD.pkl', f'{copy_to}/Gene4PD.pkl')
# copyfile('CommonData/GWASdb_SNP.pkl', f'{copy_to}/GWASdb_SNP.pkl') 
# copyfile('CommonData/LitValid_AllGenes.pkl', f'{copy_to}/LitValid_AllGenes.pkl')
# copyfile('step6_InitPreds/input/PDgenes_PDmap_noUK.pkl', f'{copy_to}/PDgenes_PDmap_noUK.pkl') 
# copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info')          

copy_from = 'step1_NMTF/output'
copy_to = f'{work_dir}/input/Genelists'                
files =  os.listdir(copy_from)    
for file in files:
    for file in files:
        if 'Geneslist_' in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

# copyfile('Data/PDmap/output/PDgenes_PDmap_noUK.pkl', f'{work_dir}/input/PDgenes_PDmap_noUK.pkl') 
# copyfile('Data/PDmap/output/PDgenes_PDmap.pkl', f'{work_dir}/input/PDgenes_PDmap.pkl') 
# copyfile('Data/PDmap/output/PDmap_PINK1.pkl', f'{work_dir}/input/PDmap_PINK1.pkl') 
copyfile('Data/DisGeNet/output/PDgenes_DGN_ALL.pkl', f'{work_dir}/input/PDgenes_DGN_ALL.pkl') 



os.chdir(work_dir) 
