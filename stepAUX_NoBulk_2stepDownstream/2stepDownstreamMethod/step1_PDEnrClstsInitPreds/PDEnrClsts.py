# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:58:57 2023

@author: kmihajlo
"""

import os, pickle
from scipy.stats import hypergeom
import numpy as np

# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p  
    

def PDenrichedClusts(G1_clust, PD_genes):
    genes = [str(item) for sublist in G1_clust for item in sublist]
    # genes = [Entrez2Sym[gene] for gene in genes]                    
    Possible_PDgenes = list(set(genes) & set(PD_genes))
                          
    Enriched_clusts = []
    PDpreds_clusts = []
    pvals = []
    Enriched_clusts_i = []
    PDgenes_cc = []
    fold_clust = []
    for i in range(len(G1_clust)): 
        genes_clust = G1_clust[i]
        # genes_clust = [Entrez2Sym[str(x)] for x in genes_clust]
        PDgenes_clust = [x for x in Possible_PDgenes if x in genes_clust]                            
        
        M = len(genes)
        K = len(Possible_PDgenes)
        N = len(genes_clust)
        X = len(PDgenes_clust)
        try:
            fold = (X/N)/(K/M)
        except ZeroDivisionError:
            fold = 0
        if fold >= 1:
            # print(fold)
            pval = hypergeom.sf(X-1, M, K, N)
            if pval <= 0.05: 
                Enriched_clusts_i.append(i)
                Enriched_clusts.append(genes_clust)
                fold_clust.append(fold)
                PDpreds = [x for x in genes_clust if x not in Possible_PDgenes]
                Clst_PDgenes = [x for x in genes_clust if x in Possible_PDgenes]
                PDgenes_cc.append(Clst_PDgenes)
                
                PDpreds_clusts.append(PDpreds)
                pvals.append(pval)

    # Benjamini-Hochberg p-value correction   
    pvals_adj =  p_adjust_bh(pvals)
    indexes = []
    for i in range(len(pvals_adj)):
        if pvals_adj[i] > 0.05:
            indexes.append(i)
    
    if len(indexes) >= 1:
        for index in sorted(indexes, reverse=True):
            del PDpreds_clusts[index]
            del Enriched_clusts[index]
            del Enriched_clusts_i[index]
            del PDgenes_cc[index]
            del fold_clust[index]
                
    # perc of clusters enriched in PD genes and percentage of total PD genes these clusters capture and percentage of preds, compared to all other genes in the condition
    enr_cluster = len(Enriched_clusts)
    total_cluster = len(G1_clust)
    perc_enr_cluster = 100.*enr_cluster/total_cluster    
           
    PDgenes_cc = [item for sublist in PDgenes_cc for item in sublist]
    percPDgenes = len(PDgenes_cc)/len(Possible_PDgenes)*100
    
    fold_clust_avg = np.mean(fold_clust)
    return perc_enr_cluster, percPDgenes, Enriched_clusts, Enriched_clusts_i, fold_clust_avg



in_dir = 'input/Clusters' 
with open('input/PDgenes_DGN_ALL.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)


 
# for cc in os.listdir(in_dir):
#     for file in os.listdir(f'{in_dir}/{cc}'):       
#         save_dir = f'output/{cc}'
#         if not os.path.exists(save_dir):
#             os.makedirs(save_dir) 

#         with open(f'{root}/{file}', 'rb') as handle:
#             G1_clust = pickle.load(handle)
           
#         perc_enr_cluster, percPDgenes, Enriched_clusts, Enriched_clusts_i,_ = PDenrichedClusts(G1_clust, PD_genes)
#         save_file1 = f'{save_dir}/{cc}_PDEnrichClusts.pkl'
#         save_file2 = f'{save_dir}/{cc}_PDEnrichClusts_indices.pkl'
#         with open(save_file1, 'wb') as handle:
#             pickle.dump(Enriched_clusts, handle)   
#         with open(save_file2, 'wb') as handle:
#             pickle.dump(Enriched_clusts_i, handle)                  



runs = 10
net_comb = 'P_1.0_C_1.0_G_0.0_M_10.0_E_10.0'


for cc in os.listdir(in_dir):
    print(cc)
    percPDgenes_cc = []
    perc_enr_cluster_cc = []
    
    for file in os.listdir(f'{in_dir}/{cc}'):       
        run = file.split('_')[0]
        print(run)            
        with open(f'{in_dir}/{cc}/{file}', 'rb') as handle:
            G1_clust = pickle.load(handle)                
        perc_enr_cluster, percPDgenes,_ ,_,_ = PDenrichedClusts(G1_clust, PD_genes)
                       
        percPDgenes_cc.append(percPDgenes)
        perc_enr_cluster_cc.append(perc_enr_cluster)
           
    tmp = max(percPDgenes_cc)
    index_max = percPDgenes_cc.index(tmp)
    print(tmp, index_max, perc_enr_cluster_cc[index_max])

    save_dir = f'output/{cc}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    clst_file = f'{index_max}_G1_clsts_{net_comb}.pkl'
    
    with open(f'input/Clusters/{cc}/{clst_file}', 'rb') as handle:
        G1_clust = pickle.load(handle)
    perc_enr_cluster, percPDgenes, Enriched_clusts, Enriched_clusts_i,_ = PDenrichedClusts(G1_clust, PD_genes)


    save_file1 = f'{save_dir}/{index_max}_PDEnrichClusts.pkl'
    save_file2 = f'{save_dir}/{index_max}_PDEnrichClusts_indices.pkl'

    with open(save_file1, 'wb') as handle:
        pickle.dump(Enriched_clusts, handle)   
    with open(save_file2, 'wb') as handle:
        pickle.dump(Enriched_clusts_i, handle)  


            
            
            
            
            
            
            
            
            
            