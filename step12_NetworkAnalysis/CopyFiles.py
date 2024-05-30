# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 14:12:56 2023

@author: kmihajlo
"""

from shutil import copyfile
import os,sys

wd = str(sys.argv[1])

if not os.path.exists(f'{wd}/input'):
    os.makedirs(f'{wd}/input') 
if not os.path.exists(f'{wd}/output'):
    os.makedirs(f'{wd}/output')
    

copy_from = 'step7_TGM_Predictions/output'
copy_to = f'{wd}/input'  

for days_set in os.listdir(copy_from):   
    print(days_set)
    for perc_case in os.listdir(f'{copy_from}/{days_set}'):         
        if os.path.isdir(f'{copy_from}/{days_set}/{perc_case}'):
            print(perc_case)              
            copy_to = f'{wd}/input/Predictions/{days_set}/{perc_case}'  
            if not os.path.exists(copy_to):
                 os.makedirs(copy_to)    
            for file in os.listdir(f'{copy_from}/{days_set}/{perc_case}'):
                if 'GM_Length' in file and file.endswith('.csv'):
                    print(file)
                    copyfile(f'{copy_from}/{days_set}/{perc_case}/{file}', f'{copy_to}/{file}')
                    
copy_to = f'{wd}/input'  

copyfile('step11_ComparisonDEGsDAPs/output/D8_D18_D25_D32_D37/DAPs_D8_D18_D25_D32_D37_1FC.csv', f'{copy_to}/DAPs_D8_D18_D25_D32_D37_1FC.csv')
copyfile('step11_ComparisonDEGsDAPs/output/D8_D18_D25_D32_D37/DEGs_D8_D18_D25_D32_D37_0.5FC.csv', f'{copy_to}/DEGs_D8_D18_D25_D32_D37_0.5FC.csv')



copy_from = 'step1_CreateNetworks/input/MolecularNetworks'
copy_to = f'{wd}/input/MolecularNetworks' 
if not os.path.exists(copy_to):
    os.makedirs(copy_to) 
        
copyfile(f'{copy_from}/Human_PPI_General.edgelist', f'{copy_to}/Human_PPI_General.edgelist')
copyfile('step1_CreateNetworks/output/Entrez2Sym.pkl', f'{wd}/input/Entrez2Sym.pkl')



copy_from = 'step2_NMTF_WeightedNets/output'
copy_to = f'{wd}/input/Genelists'                
files =  os.listdir(copy_from)    
for file in files:
    for file in files:
        if 'Geneslist_' in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')

copy_from = 'Data/ParsingGeneAnnotDBs/output'
copy_to = f'{wd}/input'     
copyfile(f'{copy_from}/HSA_GO-BP.lst', f'{copy_to}/HSA_GO-BP.lst')
copyfile(f'{copy_from}/HSA_Reactome_Pathways.lst', f'{copy_to}/HSA_Reactome_Pathways.lst')
copyfile(f'{copy_from}/HSA_Reactome_Pathways_meaning.lst', f'{copy_to}/HSA_Reactome_Pathways_meaning.lst')
copyfile(f'{copy_from}/Reactome_Group.csv', f'{copy_to}/Reactome_Group.csv')


copyfile('Data/ParsingGeneAnnotDBs/input/go-basic.obo', f'{copy_to}/go-basic.obo')
copyfile('step1_CreateNetworks/output/Entrez2Sym.pkl', f'{copy_to}/Entrez2Sym.pkl')



os.chdir(wd) 
