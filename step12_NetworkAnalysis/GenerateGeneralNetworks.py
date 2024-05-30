# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 15:30:40 2023

@author: Katarina
"""
import networkx as nx
import pickle, os, sys
import pandas as pd

wd = str(sys.argv[1])

with open(f'{wd}/input/Entrez2Sym.pkl','rb') as f:
    Entrez2Sym = pickle.load(f)   

all_genes = []
for file in os.listdir(f'{wd}/input/Genelists'):
    genes = pd.read_csv(f'{wd}/input/Genelists/{file}', header = 0)
    genes = genes['genes'].to_list()
    all_genes.append(genes)
common_genes = set.intersection(*map(set,all_genes))
union_genes =  set().union(*all_genes)
with open(f'{wd}/output/All_genes_union.txt', 'w') as f:
    for gene in union_genes:
        f.write(f'{gene}\n')


gen_nets_dir = f'{wd}/input/MolecularNetworks'

PPI = nx.read_edgelist(f'{gen_nets_dir}/Human_PPI_General.edgelist')
PPI = nx.relabel_nodes(PPI, Entrez2Sym)
nx.write_edgelist(PPI, f'{wd}/output/PPI_General_Biogrid_GeneSym.edgelist', data=False)


