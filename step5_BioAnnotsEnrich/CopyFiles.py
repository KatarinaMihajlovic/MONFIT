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
copyfile(f'{copy_from}/HSA_Kegg_Pathways.lst', f'{copy_to}/HSA_Kegg_Pathways.lst')
copyfile(f'{copy_from}/HSA_Reactome_Pathways.lst', f'{copy_to}/HSA_Reactome_Pathways.lst')
copyfile('Data/ParsingGeneAnnotDBs/input/go-basic.obo', f'{copy_to}/go-basic.obo')
copyfile(f'{copy_from}/HSA_KEGG_Pathways_meaning.lst', f'{copy_to}/HSA_KEGG_Pathways_meaning.lst')
copyfile(f'{copy_from}/HSA_Reactome_Pathways_meaning.lst', f'{copy_to}/HSA_Reactome_Pathways_meaning.lst')

copyfile('step1_CreateNetworks/output/Entrez2Sym.pkl', f'{copy_to}/Entrez2Sym.pkl')

