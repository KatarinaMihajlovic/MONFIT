# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 16:57:53 2023

@author: kmihajlo
"""

from shutil import copyfile
import os,sys

wd = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

if not os.path.exists(f'{wd}/input'):
    os.makedirs(f'{wd}/input') 
if not os.path.exists(f'{wd}/output'):
    os.makedirs(f'{wd}/output')
    
copy_from = 'step4_TGM_Predictions/output'
copy_to = f'{wd}/input'  

for days_set in os.listdir(copy_from):   
    for perc_case in os.listdir(f'{copy_from}/{days_set}'):
        if '.jpg' not in perc_case and '.csv' not in perc_case:
            copy_to = f'{wd}/input/Predictions/{days_set}/{perc_case}'  
            if not os.path.exists(copy_to):
                 os.makedirs(copy_to)                     
            for file in os.listdir(f'{copy_from}/{days_set}/{perc_case}'):
                if '.csv' in file and 'GM_Length' in file:
                    copyfile(f'{copy_from}/{days_set}/{perc_case}/{file}', f'{copy_to}/{file}')
           


copy_from = 'step1_NMTF/output'
copy_to = f'{wd}/input/Genelists'                
files =  os.listdir(copy_from)    
for file in files:
    if 'Geneslist_' in file:
        if not os.path.exists(copy_to):
            os.makedirs(copy_to)
        copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


copyfile('Data/LitValid_AllGenes/output/LitValid_AllGenes.pkl', f'{wd}/input/LitValid_AllGenes.pkl') 
#copyfile('Data/PDmap/output/PDgenes_PDmap.pkl', f'{wd}/input/PDgenes_PDmap.pkl') 
copyfile('Data/DisGeNet/output/PDgenes_DGN_ALL.pkl', f'{wd}/input/PDgenes_DGN_ALL.pkl') 
copyfile('Data/Parse_Gene4PD/output/Gene4PD.pkl', f'{wd}/input/Gene4PD.pkl') 

copyfile('Data/Homo_sapiens.gene_info', f'{wd}/input/Homo_sapiens.gene_info') 


os.chdir(wd) 

