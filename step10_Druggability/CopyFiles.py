# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:45:26 2023

@author: bscuser
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


        
copy_from = 'step2_NMTF_WeightedNets/output'
copy_to = f'{wd}/input/Genelists'                
files =  os.listdir(copy_from)    
for file in files:
    for file in files:
        if 'Geneslist_' in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')



copyfile('Data/ParsingDrugbank/output/DTI_invest_approv.csv', f'{wd}/input/DTI_invest_approv.csv') 
