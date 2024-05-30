# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 12:09:44 2023

@author: kmihajlo
"""

import pandas as pd
import os


enzyme_x = []
substrates_x = []
products_x = []
genes_id_x = []
# genes_name_x = []

human_enzymes = os.listdir('input/human_enzymes')
for enzyme in human_enzymes:
    with open(f'input/human_enzymes/{enzyme}') as f:
        lines = f.readlines()

    substrates = []
    products = []
    genes_id = []
    
    i_prod = [i for i, s in enumerate(lines) if 'PRODUCT' in s][0]
    i_cpd = [i for i, s in enumerate(lines) if 'CPD:' in s]
    i_hsa = [i for i, s in enumerate(lines) if 'HSA:' in s][0]    
    
    # substrates and products
    for i in i_cpd:
       cpd = lines[i].split(':')[1].split(']')[0]    
       if i < i_prod:
           substrates.append(cpd)
       else:
           products.append(cpd)
           
    genes_lines = lines[i_hsa].strip().split(' ')
    genes_lines = [gene for gene in genes_lines if gene != 'GENES' and gene != '' and gene != 'HSA:']    
    genes_ids = [gene.split('(')[0] for gene in genes_lines]
    
    for j in range(len(genes_ids)):
        genes_id.append(genes_ids[j])
        
    enzyme_x.append(enzyme)
    substrates_x.append(substrates)
    products_x.append(products)
    genes_id_x.append(genes_id)

if not os.path.exists('output'):
    os.makedirs('output')         

enzyme_subprodgene = pd.DataFrame(list(zip(enzyme_x,substrates_x,products_x,genes_id_x)),columns=['enzyme','substrates','products','genes_id'])     
enzyme_subprodgene = enzyme_subprodgene.set_index('enzyme')        
enzyme_subprodgene.to_csv('output/enzyme_subprodgene_KEGG.csv',index=True)   



        
        
        
        
        
        
        
        
        
        
        