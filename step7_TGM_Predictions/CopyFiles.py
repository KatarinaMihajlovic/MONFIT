# -*- coding: utf-8 -*-
"""
Created on Thu Oct  14 12:48:56 2021

@author: kmihajlo
"""

from shutil import copyfile
import os,sys

wd = str(sys.argv[1])

if not os.path.exists(f'{wd}/input'):
    os.makedirs(f'{wd}/input') 
if not os.path.exists(f'{wd}/output'):
    os.makedirs(f'{wd}/output')

copy_from = 'step2_NMTF_WeightedNets/output'
copy_to = f'{wd}/input/Genelists'                
files =  os.listdir(copy_from)    
for file in files:
    for file in files:
        if 'Geneslist_' in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')


copyfile('Data/DisGeNet/output/PDgenes_DGN_ALL.pkl', f'{wd}/input/PDgenes_DGN_ALL.pkl') 


