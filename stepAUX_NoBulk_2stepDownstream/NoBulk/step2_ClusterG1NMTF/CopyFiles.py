# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""

from shutil import copyfile
import os

work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 



if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')



copy_from = 'step1_NMTF/output'
files =  os.listdir(copy_from)    
for file in files:
    copy_to = f'{work_dir}/input/Genelists'                
    for file in files:
        if 'Geneslist_' in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')
 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
    
copyfile('Data/ParsingGeneAnnotDBs/input/go-basic.obo', f'{work_dir}/input/go-basic.obo')
copyfile('Data/PDmap/output/PDmap_BasePathways_noSBUK.lst', f'{work_dir}/input/PDmapPaths_noSBUK.lst') 

os.chdir(work_dir) 
