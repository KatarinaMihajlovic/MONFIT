# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 11:01:11 2023

@author: kmihajlo
"""

# Python program to read
# json file

import pickle
import json
    
with open('input/pd_map_HGNC_annotations.json', 'r', errors='ignore') as read_file:
    pd_map_HGNC_annotations = json.load(read_file)

# base pathways
PDmap_BasePathways = {}
for path in pd_map_HGNC_annotations:
    genes_path = set()
    for i in range(len(pd_map_HGNC_annotations[path])):
        gene = pd_map_HGNC_annotations[path][i]['resource']
        genes_path.add(gene)
    PDmap_BasePathways[path] = genes_path 
with open('output/PDmap_BasePathways.pkl', 'wb') as handle:
    pickle.dump(PDmap_BasePathways, handle)  


PDmap_BasePathways_output = open('output/PDmap_BasePathways.lst','w')
for path in PDmap_BasePathways:
    for gene in PDmap_BasePathways[path]:
        PDmap_BasePathways_output.write(gene + '\t' + path + '\n')
PDmap_BasePathways_output.close()


PDmap_BasePathways.pop('Parkinsons UK Gene Ontology genes')
PDmap_BasePathways.pop('Scrapbook')

PDmap_BasePathways_output = open('output/PDmap_BasePathways_noSBUK.lst','w')
for path in PDmap_BasePathways:
    for gene in PDmap_BasePathways[path]:
        PDmap_BasePathways_output.write(gene + '\t' + path + '\n')
PDmap_BasePathways_output.close()

with open('output/PDmap_BasePathways_noSBUK.pkl', 'wb') as handle:
    pickle.dump(PDmap_BasePathways, handle)  









