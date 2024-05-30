# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""

from shutil import copyfile
import os, sys

wd = str(sys.argv[1])

if not os.path.exists(f'{wd}/input'):
    os.makedirs(f'{wd}/input') 
if not os.path.exists(f'{wd}/output'):
    os.makedirs(f'{wd}/output')



copy_from = 'step2_NMTF_WeightedNets/output'
files =  os.listdir(copy_from)    
for file in files:
    copy_to = f'{wd}/input/Genelists'                
    for file in files:
        if 'Geneslist_' in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')
    
copyfile('Data/PDmap/output/PDmap_BasePathways_noSBUK.lst', f'{wd}/input/PDmapPaths_noSBUK.lst') 

