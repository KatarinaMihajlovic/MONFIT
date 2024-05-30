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

Metabolomics_merged = pd.read_csv('input/Metabolomics/Metabolomics_merged_KEGGid.csv', index_col=0)
# print(len(Metabolomics_merged.columns))
# plot_hist(Metabolomics_merged.values.flatten(), title='x', logs = False)
enzyme_subprodgene_KEGG = pd.read_csv('input/MolecularNetworks/enzyme_subprodgene_KEGG.csv', index_col=0)

Proteomics = pd.read_csv('input/Proteomics/Proteomicslog.csv', index_col=0)
# plot_hist(Proteomics.values.flatten(), title='x', logs = False)

metab_median_df = Metabolomics_merged.median(axis=1)
prot_min_df = Proteomics[Proteomics > 0].min(axis=1)

# mappingMets = pd.read_excel('input/Metabolomics/MappingMetabolite names_Katarina_SkupinLab.xlsx')
# mapDict = dict(zip(mappingMets.KEGG,mappingMets.Query))
with open('input/Metabolomics/KEGGID2metab.pkl','rb') as f:
    mapDict = pickle.load(f)    

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


def MI_enzymeCentr_w(genes_cc, cc, Metabolomics_merged, enzyme_subprodgene_KEGG, metab_median_df, weighted = False):
    metab_median = metab_median_df.loc[cc]
    
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
    
    # if weighted == True:
    metabolites = list(Metabolomics_merged.columns)   
    # else:
    #     subs_flat = list(set([item for sublist in substrates_x for item in sublist]))
    #     prods_flat = list(set([item for sublist in products_x for item in sublist]))
    #     metabolites = list(set(subs_flat) & set(prods_flat))
    uniquemets = set()

    MI_net_metab = nx.Graph()
    for prods_id in range(len(products_x)):
        for prod in products_x[prods_id]:
            for subs_id in range(len(substrates_x)):
                # print(substrates_x[subs_id])
                if prod in substrates_x[subs_id]:
                    metabolite = prod
                    genes_prods = genes_id_x[prods_id] 
                    genes_subs = genes_id_x[subs_id]           
                    for gene_prod in genes_prods:
                        gene_prod = Entrez2Sym[gene_prod]
                        for gene_sub in genes_subs:
                            gene_sub = Entrez2Sym[gene_sub]
                            if gene_prod in genes_cc and gene_sub in genes_cc and gene_prod != gene_sub:  
                                if weighted != True:
                                    MI_net_metab.add_edge(gene_prod, gene_sub, weight=1)
                                else:                               
                                    if MI_net_metab.has_edge(gene_prod, gene_sub) != True:
                                        if prod in metabolites:
                                            uniquemets.add(metabolite)
                                            weight = Metabolomics_merged.loc[cc,metabolite]
                                            MI_net_metab.add_edge(gene_prod, gene_sub, weight=weight)   
                                        else:
                                            MI_net_metab.add_edge(gene_prod, gene_sub, weight=metab_median)   
                                    else: 
                                        w = MI_net_metab.get_edge_data(gene_prod, gene_sub)
                                        if prod in metabolites:
                                            uniquemets.add(metabolite)
                                            weight = Metabolomics_merged.loc[cc,metabolite]
                                            weight = w['weight'] + weight
                                            MI_net_metab.add_edge(gene_prod, gene_sub, weight=weight)   
                                        else:
                                            weight = w['weight'] + metab_median
                                            MI_net_metab.add_edge(gene_prod, gene_sub, weight=weight)   
                                            
                                        
    # if weighted == True:    
    #     MI_df = nx.to_pandas_adjacency(MI_net_metab) 
    #     MI_df_vals = MI_df.values.flatten()
    #     plot_hist(MI_df_vals, title = f'{cc} - MI', logs=False)
    return MI_net_metab, list(uniquemets)

def scale_net_weights(net):
    net_df = nx.to_pandas_adjacency(net)
    net_sum = np.sum(net_df.values.flatten())  
    net_df_sumdiv = net_df/net_sum
    net_sumdiv = nx.from_pandas_adjacency(net_df_sumdiv)  
    return net_sumdiv


 


    

# use transcriptomics data to filter all matrices for these expressed genes, and for genes that are in PPI
transc_data_dir = 'input/Transcriptomics'  
gen_nets_dir = 'input/MolecularNetworks'


PPI = nx.read_edgelist(f'{gen_nets_dir}/Human_PPI_General.edgelist')
PPI = nx.relabel_nodes(PPI, Entrez2Sym)
nx.write_edgelist(PPI, 'output/PPI_General.edgelist')
genes_PPI = list(PPI.nodes())
print(f'PPI \n genes: {PPI.number_of_nodes()} \n interactions: {PPI.number_of_edges()}\n density: {nx.density(PPI)}\n\n')
# print(len(genes_PPI))
        
ccs = []
All_genes = set()
metabs = []

for file in os.listdir(transc_data_dir):
    print(file)
    cc = file.split('.')[0].split('_')[1] + '_' + file.split('.')[0].split('_')[2]
    ccs.append(cc)
    sd = create_directory(f'output/{cc}')


    # read EXP
    result = pyreadr.read_r(f'{transc_data_dir}/{file}')
    EXP_cc = result[None]
    EXP_cc = EXP_cc.loc[(EXP_cc.sum(axis=1) != 0), (EXP_cc.sum(axis=0) != 0)]
    expressed_genes = EXP_cc.index.tolist()
    gene_overlap_EXP = list(set(expressed_genes) & set(genes_PPI))
    print(len(gene_overlap_EXP))
    
    ### E    
    print('E')
    EXP_cc_filter = EXP_cc.loc[gene_overlap_EXP]
 
    EXP_cc_filter = EXP_cc_filter.T
    adata = anndata.AnnData(X=EXP_cc_filter)
    adata.obs['total_counts'] = adata.X.sum(axis=1)
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)

    # log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
    axes[1].set_title("Shifted logarithm")
    plt.show()
    
    EXP_cc_normalized = pd.DataFrame(adata.layers["log1p_norm"], index = EXP_cc_filter.index, columns = EXP_cc_filter.columns)
    plot_hist(EXP_cc_normalized.values.flatten(), f'E_{cc}', sd, f'E_{cc}',logs = False)       

    # save E
    EXP_cc_normalized = EXP_cc_normalized.T
    EXP_cc_normalized.to_csv(f'{sd}/E_{cc}_normalized.csv' , header = True, index = True)        
    
    cc_genes = EXP_cc_normalized.index.tolist()
    All_genes.update(cc_genes)
    
    removed_genes = list(set(EXP_cc.index.to_list()) - set(gene_overlap_EXP))
    with open(f'{sd}/EXP_{cc}_RemovedGenes.txt', 'w') as f:
        for rg in removed_genes:
            f.write(f"{rg}\n")
       
   
    ### Molecular Networks
    Nets_info = open(f'{sd}/Network_Statistic_{cc}.txt','w')

    ### PPI
    print('PPI')
    # create a condition specific PPI graph using genes measured by  E
    PPI_EXP = nx.Graph(PPI.subgraph(gene_overlap_EXP))
    PPI_EXP.number_of_edges()
    PPI_EXP_df = nx.to_pandas_adjacency(PPI_EXP)          

    # weight all edges with product of min values from proteomics
    prot_min = prot_min_df.loc[cc]
    PPI_dfmin = PPI_EXP_df.replace(1,prot_min*prot_min)
    # nx.write_edgelist(nx.from_pandas_adjacency(PPI_dfmin), f'{sd}/PPI_dfmin_{cc}.edgelist', data=['weight'])  

    # create a PPI subgraph  of genes present in proteomics data
    cc_Prot = pd.DataFrame(Proteomics.loc[cc], columns=[cc])  
    gene_overlap_prot = list(set(gene_overlap_EXP) & set (cc_Prot.index))
    gene_overlap_prot = sorted(gene_overlap_prot)
    cc_Prot = cc_Prot.loc[gene_overlap_prot]               
    PPI_prot = nx.Graph(PPI.subgraph(gene_overlap_prot))
    PPI_prot.number_of_edges()

    PPIprot_df = nx.to_pandas_adjacency(PPI_prot)
    PPIprot_df = PPIprot_df[gene_overlap_prot]
    PPIprot_df = PPIprot_df.reindex(gene_overlap_prot)
   
    # weight edges between measured proteins with product of measured protomics abundance between the two proteins        
    p = cc_Prot.values
    ppT= np.matmul(p,p.T)  
    PPIprot_vals = PPIprot_df.values
    PPI_ppT = PPIprot_vals * ppT  

    PPI_ppT_df = pd.DataFrame(PPI_ppT, columns=gene_overlap_prot, index=gene_overlap_prot) 
    # nx.write_edgelist(nx.from_pandas_adjacency(PPI_ppT_df), f'{sd}/PPI_ppT_df_{cc}.edgelist', data=['weight'])  
    # PPI_ppT_df.to_csv( f'{sd}/PPI_ppT_df_{cc}.csv',index=False)


    PPI_ppT_allexp_df = PPI_ppT_df.combine_first(PPI_dfmin)
    PPI_ppT_allexp_net = nx.from_pandas_adjacency(PPI_ppT_allexp_df)       
    # PPI_ppT_allexp_net.number_of_edges()

    # PPI_ppT = PPI_ppT_allexp_df.values
    # plot_hist(PPI_ppT.flatten(), f'{cc} - PPI*ppT', sd, f'PPI_ppT_{cc}', logs=False)
            
    nx.write_edgelist(PPI_ppT_allexp_net, f'{sd}/PPI_ppTW_{cc}.edgelist', data=['weight'])           
    Nets_info.write(f'PPI \n genes: {PPI_ppT_allexp_net.number_of_nodes()} \n interactions: {PPI_ppT_allexp_net.number_of_edges()}\n density: {nx.density(PPI_ppT_allexp_net)}\n\n')
  


    ### COEX
    net = 'COEX'
    COEX = CreateNetwork(net, gene_overlap_EXP, gen_nets_dir)    
    COEX_df = nx.to_pandas_adjacency(COEX)
    nx.write_edgelist(COEX, f'{sd}/{net}_{cc}.edgelist') 
    
    # plot_hist(COEX_df.values.flatten(), f'{cc} - COEX', sd, f'COEX_{cc}', logs = False)

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
    MI, uniquemets = MI_enzymeCentr_w(gene_overlap_EXP, cc, Metabolomics_merged, enzyme_subprodgene_KEGG, metab_median_df, weighted = True)
    MI_df = nx.to_pandas_adjacency(MI)
    net_name =  f'{net}_MetAbundWeight_{cc}.edgelist'
    nx.write_edgelist(MI, f'{sd}/{net_name}', data=['weight'])      
    # plot_hist(MI_df.values.flatten(), f'{cc} - MI', sd,  f'MI_{cc}', logs = False)  
    metabs.append(uniquemets)        
    Nets_info.write(f'{net} \n genes: {MI.number_of_nodes()} \n interactions: {MI.number_of_edges()}\n density: {nx.density(MI)}\n\n')        

    Nets_info.close()   

        
                
All_gene_df = pd.DataFrame(All_genes, columns = ['genes'])
All_gene_df.to_csv('output/All_genes.csv',index=False)
                




