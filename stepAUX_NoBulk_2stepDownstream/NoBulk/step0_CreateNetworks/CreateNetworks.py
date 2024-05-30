# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:10:47 2023

@author: kmihajlo
"""

import os,pyreadr
import pandas as pd
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import anndata
import scanpy as sc
import seaborn as sns

def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    Entrez2Sym_synonyms = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        organism = lspt[0]
        if organism == '9606':
            Symbol = lspt[2]
            synonyms = lspt[4]
            if synonyms != '-':    
                if '|' in synonyms:
                    synonyms = synonyms.split('|')
                    Entrez2Sym_synonyms[lspt[1]] = synonyms
                else:
                    Entrez2Sym_synonyms[lspt[1]] = synonyms
    
            Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym,Entrez2Sym_synonyms)

Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym,Entrez2Sym_synonyms = Entrez2Symbols(Entrez2Symbol_file)   
with open('output/Entrez2Sym.pkl', 'wb') as f:
    pickle.dump(Entrez2Sym, f)

enzyme_subprodgene_KEGG = pd.read_csv('input/MolecularNetworks/enzyme_subprodgene_KEGG.csv', index_col=0)

def create_directory(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name) 
    return dir_name


def CreateNetwork(GeneralNet, gene_overlap, gen_nets_dir, Entrez2Sym = Entrez2Sym):
    print(GeneralNet)
    gennet = nx.read_edgelist(f'{gen_nets_dir}/Human_{GeneralNet}_General.edgelist')

    gennet = nx.relabel_nodes(gennet, Entrez2Sym)
    gennet = gennet.subgraph(gene_overlap)
    return gennet

def plot_hist(a, title, save_dir, figname,logs = True):
    a = a[a != 0]
    # a = np.sort(a)[::-1]
    fig, ax = plt.subplots(figsize=(8, 5))
    plt.hist(a, bins = 100) 
    # plt.axvline(x = v_line, color = 'r')
    if logs == True:
        plt.yscale('log')
    plt.title(title, fontsize=22) 
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)

    plt.savefig(f'{save_dir}/{figname}_hist.jpg', dpi = 350, format='jpg')  
    plt.show()
    plt.close()

def MI_EnzymeCentr(enzyme_subprodgene_KEGG, genes_cc):
    substrates_x = []
    substrates_adj = enzyme_subprodgene_KEGG['substrates'].to_list()
    for sub in substrates_adj:
        subs = sub.split("'")
        subs = [x for x in subs if 'C' in x]
        substrates_x.append(subs)
    
    products_x = []    
    products_adj = enzyme_subprodgene_KEGG['products'].to_list()
    for prod in products_adj:
        prods = prod.split("'")
        prods = [x for x in prods if 'C' in x]
        products_x.append(prods)
        
    genes_id_x = []
    genes_id_adj = enzyme_subprodgene_KEGG['genes_id'].to_list()
    for genes in genes_id_adj:
        genes = genes.split("'")
        genes = [x for x in genes if x != ', ' and x!=']' and x!='[']
        genes_id_x.append(genes)
    
    MI_net = nx.Graph()
    for prods_id in range(len(products_x)):
        for prod in products_x[prods_id]:
            for subs_id in range(len(substrates_x)):
                if prod in substrates_x[subs_id]:
                    genes_prod = genes_id_x[prods_id]
                    genes_subs = genes_id_x[subs_id]
                    for gene_prod in genes_prod:
                        gene_prod = Entrez2Sym[gene_prod]
                        for gene_sub in genes_subs:
                            gene_sub = Entrez2Sym[gene_sub]
                            if gene_prod in genes_cc and gene_sub in genes_cc and gene_prod != gene_sub:  
                                MI_net.add_edge(gene_prod, gene_sub)
    return MI_net
    
   

# def MI_enzymeCentr_w(genes_cc, cc, Metabolomics_merged, enzyme_subprodgene_KEGG, metab_median_df, weighted = False):
#     metab_median = metab_median_df.loc[cc]
    
#     substrates_x = []
#     substrates_adj = enzyme_subprodgene_KEGG['substrates'].to_list()
#     for sub in substrates_adj:
#         subs = sub.split("'")
#         subs = [x for x in subs if 'C' in x]
#         substrates_x.append(subs)
    
#     products_x = []    
#     products_adj = enzyme_subprodgene_KEGG['products'].to_list()
#     for prod in products_adj:
#         prods = prod.split("'")
#         prods = [x for x in prods if 'C' in x]
#         products_x.append(prods)
        
#     genes_id_x = []
#     genes_id_adj = enzyme_subprodgene_KEGG['genes_id'].to_list()
#     for genes in genes_id_adj:
#         genes = genes.split("'")
#         genes = [x for x in genes if x != ', ' and x!=']' and x!='[']
#         genes_id_x.append(genes)
    
#     # if weighted == True:
#     metabolites = list(Metabolomics_merged.columns)   
#     # else:
#     #     subs_flat = list(set([item for sublist in substrates_x for item in sublist]))
#     #     prods_flat = list(set([item for sublist in products_x for item in sublist]))
#     #     metabolites = list(set(subs_flat) & set(prods_flat))
#     uniquemets = set()

#     MI_net_metab = nx.Graph()
#     for prods_id in range(len(products_x)):
#         for prod in products_x[prods_id]:
#             for subs_id in range(len(substrates_x)):
#                 # print(substrates_x[subs_id])
#                 if prod in substrates_x[subs_id]:
#                     metabolite = prod
#                     genes_prods = genes_id_x[prods_id] 
#                     genes_subs = genes_id_x[subs_id]           
#                     for gene_prod in genes_prods:
#                         gene_prod = Entrez2Sym[gene_prod]
#                         for gene_sub in genes_subs:
#                             gene_sub = Entrez2Sym[gene_sub]
#                             if gene_prod in genes_cc and gene_sub in genes_cc and gene_prod != gene_sub:  
#                                 if weighted != True:
#                                     MI_net_metab.add_edge(gene_prod, gene_sub, weight=1)
#                                 else:                               
#                                     if MI_net_metab.has_edge(gene_prod, gene_sub) != True:
#                                         if prod in metabolites:
#                                             uniquemets.add(metabolite)
#                                             weight = Metabolomics_merged.loc[cc,metabolite]
#                                             MI_net_metab.add_edge(gene_prod, gene_sub, weight=weight)   
#                                         else:
#                                             MI_net_metab.add_edge(gene_prod, gene_sub, weight=metab_median)   
#                                     else: 
#                                         w = MI_net_metab.get_edge_data(gene_prod, gene_sub)
#                                         if prod in metabolites:
#                                             uniquemets.add(metabolite)
#                                             weight = Metabolomics_merged.loc[cc,metabolite]
#                                             weight = w['weight'] + weight
#                                             MI_net_metab.add_edge(gene_prod, gene_sub, weight=weight)   
#                                         else:
#                                             weight = w['weight'] + metab_median
#                                             MI_net_metab.add_edge(gene_prod, gene_sub, weight=weight)   
                                            
                                        
#     # if weighted == True:    
#     #     MI_df = nx.to_pandas_adjacency(MI_net_metab) 
#     #     MI_df_vals = MI_df.values.flatten()
#     #     plot_hist(MI_df_vals, title = f'{cc} - MI', logs=False)
#     return MI_net_metab, list(uniquemets)


def scale_net_weights(net):
    net_df = nx.to_pandas_adjacency(net)
    net_sum = np.sum(net_df.values.flatten())  
    net_df_sumdiv = net_df/net_sum
    net_sumdiv = nx.from_pandas_adjacency(net_df_sumdiv)  
    return net_sumdiv


 


    

# use transcriptomics data to filter all matrices for these expressed genes, and for genes that are in PPI
# transc_data_dir = 'input/Transcriptomics'  
gen_nets_dir = 'input/MolecularNetworks'


PPI = nx.read_edgelist(f'{gen_nets_dir}/Human_PPI_General.edgelist')
PPI = nx.relabel_nodes(PPI, Entrez2Sym)
nx.write_edgelist(PPI, 'output/PPI_General.edgelist')
genes_PPI = list(PPI.nodes())
print(f'PPI \n genes: {PPI.number_of_nodes()} \n interactions: {PPI.number_of_edges()}\n density: {nx.density(PPI)}\n\n')
# print(len(genes_PPI))
        
All_genes = set()

cell_conds = ['ND_18', 'WT_18', 'WT_32', 'ND_32', 'WT_8', 'ND_8', 'ND_25', 'WT_25', 'WT_37', 'ND_37']

for cc in cell_conds:
    print(cc)
    sd = create_directory(f'output/{cc}')
    EXP = pd.read_csv(f'input/{cc}/E_{cc}_normalized.csv', index_col=0)
    expressed_genes = EXP.index.tolist()
    gene_overlap_EXP = list(set(expressed_genes) & set(genes_PPI))


    ### Molecular Networks
    Nets_info = open(f'{sd}/Network_Statistic_{cc}.txt','w')

    ### PPI
    print('PPI')
    # create a condition specific PPI graph using genes measured by  E
    PPI_EXP = nx.Graph(PPI.subgraph(gene_overlap_EXP))
    nx.write_edgelist(PPI_EXP, f'{sd}/PPI_{cc}.edgelist')           
    Nets_info.write(f'PPI \n genes: {PPI_EXP.number_of_nodes()} \n interactions: {PPI_EXP.number_of_edges()}\n density: {nx.density(PPI_EXP)}\n\n')

  
    ### COEX
    net = 'COEX'
    COEX = CreateNetwork(net, gene_overlap_EXP, gen_nets_dir)    
    nx.write_edgelist(COEX, f'{sd}/{net}_{cc}.edgelist') 
    Nets_info.write(f'{net} \n genes: {COEX.number_of_nodes()} \n interactions: {COEX.number_of_edges()}\n density: {nx.density(COEX)}\n\n')        
    
    
    ### GI
    net = 'GI'
    GI = CreateNetwork(net, gene_overlap_EXP, gen_nets_dir)    
    GI_df = nx.to_pandas_adjacency(GI)
    nx.write_edgelist(GI, f'{sd}/{net}_{cc}.edgelist')         
    Nets_info.write(f'{net} \n genes: {GI.number_of_nodes()} \n interactions: {GI.number_of_edges()}\n density: {nx.density(GI)}\n\n')        
          
    # MI 
    net = 'MI'
    print(net)  
    MI=MI_EnzymeCentr(enzyme_subprodgene_KEGG, gene_overlap_EXP)
    
    # MI, uniquemets = MI_enzymeCentr_w(gene_overlap_EXP, cc, Metabolomics_merged, enzyme_subprodgene_KEGG, metab_median_df, weighted = True)
    # MI_df = nx.to_pandas_adjacency(MI)
    # net_name =  f'{net}_MetAbundWeight_{cc}.edgelist'
    nx.write_edgelist(MI, f'{sd}/{net}_{cc}.edgelist', data=['weight'])      
    # plot_hist(MI_df.values.flatten(), f'{cc} - MI', sd,  f'MI_{cc}', logs = False)  
    # metabs.append(uniquemets)        
    Nets_info.write(f'{net} \n genes: {MI.number_of_nodes()} \n interactions: {MI.number_of_edges()}\n density: {nx.density(MI)}\n\n')        

    Nets_info.close()   

        
                
All_gene_df = pd.DataFrame(All_genes, columns = ['genes'])
All_gene_df.to_csv('output/All_genes.csv',index=False)
                




