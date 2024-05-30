# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:07:18 2024

@author: Katarina
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os
import pandas as pd

work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')


copy_from = 'step4_Enrich_PDgenesPathways/output'
copy_to = f'{work_dir}/input'         
copyfile(f'{copy_from}/Enrichment_Rank.csv', f'{copy_to}/Enrichment_Rank.csv')

os.chdir(work_dir) 
