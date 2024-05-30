# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 15:44:45 2023

@author: kmihajlo
"""
import networkx as nx
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from random import sample
import os, sys
import pandas as pd
from statannotations.Annotator import Annotator
import itertools


def BasicStats(G, nodes):
    G_CPD = G.subgraph(nodes)
    
    nodes = nx.number_of_nodes(G_CPD)
    edges = nx.number_of_edges(G_CPD)
    density = nx.density(G_CPD)
    
    print(nodes, 'nodes')
    print(edges, 'edges')
    print(density, 'density')
    
    
    GLCC = G_CPD.subgraph(max(nx.connected_components(G_CPD), key=len)).copy()
    largest_cc = max(nx.connected_components(G_CPD), key=len)
    print(len(largest_cc), 'largest connected component')
    GLCC = G_CPD.subgraph(largest_cc).copy() 
    print('density: ', nx.classes.function.density(G_CPD))
    clustercoef = nx.average_clustering(GLCC)
    print('Clustering Coefficient (AVG): ', clustercoef) 
    return largest_cc

def ShortestPath2Gene(G, gene, TargetGenes):
    G2gene = nx.single_source_dijkstra(G, gene)
    # print(G2gene)
    G2gene_lens = G2gene[0]
    G2gene_lens.pop(gene)

    TGs_SP = {key:G2gene_lens[key] for key in G2gene_lens if key in TargetGenes}
    BG_SP = {key:G2gene_lens[key] for key in G2gene_lens if key not in TargetGenes}
    TGs_SP_dist = np.array([TGs_SP[key] for key in TGs_SP])
    BG_SP_dist = np.array([BG_SP[key] for key in BG_SP])
    return TGs_SP_dist, BG_SP_dist, TGs_SP, BG_SP

def ShortestPath2Gene_MULTI(G, gene, TargetGenes):
    G2gene = nx.single_source_dijkstra(G, gene)
    # print(G2gene)
    G2gene_lens = G2gene[0]
    G2gene_lens.pop(gene)

    results_union = set().union(*TargetGenes)

    TGs_SPs = []
    TGs_SP_dists = []
    for l in TargetGenes:
        TGs_SP = {key:G2gene_lens[key] for key in G2gene_lens if key in l}
        TGs_SPs.append(TGs_SP)
        TGs_SP_dist = np.array([TGs_SP[key] for key in TGs_SP])
        TGs_SP_dists.append(TGs_SP_dist)

    BG_SP = {key:G2gene_lens[key] for key in G2gene_lens if key not in results_union}
    BG_SP_dist = np.array([BG_SP[key] for key in BG_SP])
    return TGs_SP_dists, BG_SP_dist, TGs_SPs, BG_SP

def format_int_with_commas(x):
    """
    Formats an integer with commas as thousand separators.
    """
    return f"{x:,}"


def write_txt(genes, sd, save_file):
    with open(f'{sd}/{save_file}.txt', 'w') as f: 
        for gene in genes:                        
            f.write(f'{gene}\n')
        f.write('\n')
        
def SampWithReplac_TestSubgraph(graph, subgraph_nodes, gene, propertyn, num_permutations = 10000):
    subgraph_nodes_gene = subgraph_nodes + [gene]
    if propertyn == 'density':
        observed_property = nx.density(graph.subgraph(subgraph_nodes_gene))
    elif propertyn == 'betweenness_centrality':
        observed_property = nx.betweenness_centrality(graph.subgraph(subgraph_nodes_gene))
        # observed_property = observed_property['PINK1']
        
    print(propertyn, observed_property)
    nodes_G = list(graph.nodes())
    # print(nodes_G)
    n_samp = len(subgraph_nodes_gene)-1
    
    Successes = 0
    permutation_results = np.zeros(num_permutations)
    for i in range(num_permutations):
        n_rand_nodes = sample(nodes_G, n_samp) + [gene]
        if propertyn == 'density':
            RandomNodes_property = nx.density(graph.subgraph(n_rand_nodes))
        elif propertyn == 'betweenness_centrality':
            RandomNodes_property = nx.betweenness_centrality(graph.subgraph(n_rand_nodes))
            print(RandomNodes_property)
            # RandomNodes_property = RandomNodes_property['PINK1']
        permutation_results[i] = RandomNodes_property
        if RandomNodes_property >= observed_property:
            Successes+=1   
    p_value = (Successes + 1)/(num_permutations+1) 
    return p_value
   

def EnrichHyperheom(first_neigh, gene_set, first_neigh_BG, all_nodes):
    X = len(first_neigh)
    N = len(gene_set)
    K = len(first_neigh_BG) + len(first_neigh)
    M = len(all_nodes)
    fold = (X/N)/(K/M)
    pval = hypergeom.sf(X-1, M, K, N)
    print(f'1st Neighbours to {target_gene} in {network} significance: fold = {fold}, pval = {pval}')
    return pval, fold

def survey(results, category_names, df_annots, order, pairs):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = ['lightgreen', 'cornflowerblue', 'gold', 'lightgrey']
    threshold = 1


    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_ylabel('Gene group', fontsize=30) #dist_measure_l[dist_meas]
    ax.set_xlabel('Relative gene count (%)', fontsize=26)

    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        rects = ax.barh(labels, widths, left=starts, height=0.6,
                        label=colname, color=color)

        # r, g, b, _ = color
        text_color = 'white' if color == 'cornflowerblue' else 'black'
        labels_ratio = [round(float(v), 1) if v > threshold else "" for v in rects.datavalues]  
        print(labels_ratio)

        ax.bar_label(rects, labels = labels_ratio, label_type='center', color=text_color, fontsize = 20, weight='bold', padding=5)
    ax.legend(ncols=len(category_names), bbox_to_anchor=(-0.03, 1), loc='lower left', title="Shortest path length", fontsize=24, title_fontsize=26)
    plt.yticks(fontsize=26)
    
    
    annotator = Annotator(ax, pairs, data=df_annots, x='ShortPath', y='GeneGroup', order=order, orient='h')
    annotator.configure(test='Mann-Whitney-ls', text_format='star',fontsize=22, loc='outside')
    annotator.configure(comparisons_correction="BH", correction_format="replace")
    annotator.apply_and_annotate()
    return fig, ax
  
def firstneighSignificance(case, G, gene, TargetGenes):
    print(case)
    TGs_SP_dist, BG_SP_dist, TGs_SP, BG_SP = ShortestPath2Gene(G, gene, TargetGenes)
    # print(TGs_SP_dist)
    first_neigh = [x for x in TGs_SP_dist if x == 1]
    print(len(first_neigh))
    first_neigh_BG = [key for key in BG_SP_dist if key == 1]
    pval, fold = EnrichHyperheom(first_neigh, TargetGenes, first_neigh_BG, G.nodes())
    return pval, fold

    
'''MAIN CODE'''
     
# wd = str(sys.argv[1])

wd = os.getcwd()
with open(f'{wd}/output/All_genes_Union.txt') as f:
    All_genes_Union = f.read().splitlines() 
    
target_gene = 'PINK1'
network = 'PPI'   
print(network)
G = nx.read_edgelist(f'{wd}/output/PPI_General_Biogrid_GeneSym.edgelist', data = False)
G  = G.subgraph(All_genes_Union)
print(len(G.nodes()))

in_dir = f'{wd}/input/Predictions'
for days_set in os.listdir(in_dir): 
    days = days_set.split('_')
    days = [x[1:] for x in days]
                    
    All_genes_inters = []
    for file in os.listdir(f'{wd}/input/Genelists'):
        day = file.split('_')[2][:-4]
        if day in days:
            Geneslist = pd.read_csv(f'{wd}/input/Genelists/{file}')
            Geneslist = Geneslist['genes'].tolist()
            All_genes_inters.append(Geneslist)
    All_genes_inters = list(set.intersection(*map(set,All_genes_inters)))
    
    
    for perc_case in os.listdir(f'{in_dir}/{days_set}'):
        files = os.listdir(f'{in_dir}/{days_set}/{perc_case}')
        # files.append('GM_D18_D25_D32_D37.csv')
        for file in files:       

            
            MONFITpreds = pd.read_csv(f'{in_dir}/{days_set}/{perc_case}/{file}', index_col=0)
            MONFITpreds = list(MONFITpreds.index)
            # ###  DEGs, DAPs
            ## DEGs
            DEGs_df = pd.read_csv(f'{wd}/input/DEGs_D8_D18_D25_D32_D37_0.5FC.csv', index_col=0)
            DEGs = list(DEGs_df.index)
            DEGs = list(set(All_genes_inters) & set(DEGs))            

            ## DAPs
            DAPs_df = pd.read_csv(f'{wd}/input/DAPs_D8_D18_D25_D32_D37_1FC.csv', index_col=0)
            DAPs = list(DAPs_df.index)
            DAPs = list(set(All_genes_Union) & set(DAPs))                      

            print(len(DEGs), len(DAPs))
               
            
            sd = f'{wd}/output/{days_set}_{perc_case}'
            if not os.path.exists(sd):
                os.makedirs(sd)
                
            Statistics = open(f'{sd}/Statistics.txt','w')
            
            Statistics.write('DENSITY\n')
            prop = 'density'   
            Statistics.write(f'Is the {prop} of the subgraph of the {network} induced by a gene set with {target_gene} higher than expected at random?\n')
            Statistics.write('Test: Sampling with replacement\n')

            print(prop, 'pvalue - is density bigger than random (sampling with replacement)')
            print('MONFIT')
            subgraph_nodes_gene = MONFITpreds + [target_gene]
            dens_PDpreds = nx.density(G.subgraph(subgraph_nodes_gene))
            p_value = SampWithReplac_TestSubgraph(G, MONFITpreds, target_gene, prop)
            print(dens_PDpreds, p_value)
            Statistics.write('MONFIT\n')
            Statistics.write(f'density = {dens_PDpreds}\n')            
            Statistics.write(f'pvalue = {p_value}\n\n')

            
            print('DEGs')
            subgraph_nodes_gene = DEGs + [target_gene]
            dens_DEGs = nx.density(G.subgraph(subgraph_nodes_gene))
            p_value = SampWithReplac_TestSubgraph(G, DEGs, target_gene, prop)
            print(dens_DEGs, p_value)
            Statistics.write('DEGs\n')
            Statistics.write(f'density = {dens_DEGs}\n')            
            Statistics.write(f'pvalue = {p_value}\n\n')
            
            
            print('DAPs')
            subgraph_nodes_gene = DAPs + [target_gene]
            dens_DAPs = nx.density(G.subgraph(subgraph_nodes_gene))
            p_value = SampWithReplac_TestSubgraph(G, DAPs, target_gene, prop)
            print(dens_DAPs, p_value)
            Statistics.write('DAPs\n')
            Statistics.write(f'density = {dens_DAPs}\n')            
            Statistics.write(f'pvalue = {p_value}\n\n')
            
            

            # ######### boxplot       
            TGs_SP_dist, BG_SP_dist, TGs_SP, BG_SP = ShortestPath2Gene_MULTI(G, target_gene, [MONFITpreds, DEGs, DAPs])
            avgs_SP = [np.mean(x) for x in TGs_SP_dist]
            avg_BG = np.mean(BG_SP_dist)
            print(f'AVG shortest path (in hops) to {target_gene} in {network}')
            print('MONFIT', avgs_SP[0])
            print('DEGs', avgs_SP[1])
            print('DAPs', avgs_SP[2])
            print('BG', avg_BG)
            
            Preds_SP_cnt = np.bincount(TGs_SP_dist[0])/len(TGs_SP_dist[0])*100
            BG_SP_cnt = np.bincount(BG_SP_dist)/len(BG_SP_dist)*100
            DEGs_SP_cnt = np.bincount(TGs_SP_dist[1])/len(TGs_SP_dist[1])*100
            DAPs_SP_cnt = np.bincount(TGs_SP_dist[2])/len(TGs_SP_dist[2])*100

            Preds_SP_cnt = np.delete(Preds_SP_cnt, 0)
            DEGs_SP_cnt = np.delete(DEGs_SP_cnt, 0)
            DAPs_SP_cnt = np.delete(DAPs_SP_cnt, 0)
            BG_SP_cnt = np.delete(BG_SP_cnt, 0)
            
            Preds_SPs = np.unique(TGs_SP_dist[0])
            DEGs_SPs = np.unique(TGs_SP_dist[1])
            DAPs_SPs = np.unique(TGs_SP_dist[2])
            BG_SPs = np.unique(BG_SP_dist)
                
  
            ### stacked bar chart          
            Preds_SPs_dict = dict(zip(Preds_SPs, Preds_SP_cnt))
            DEGs_SPs_dict = dict(zip(DEGs_SPs, DEGs_SP_cnt))
            DAPs_SPs_dict = dict(zip(DAPs_SPs, DAPs_SP_cnt))
            BG_SPs_dict = dict(zip(BG_SPs, BG_SP_cnt))
            for key in BG_SPs_dict:
                if key not in Preds_SPs_dict:
                    Preds_SPs_dict[key] = 0
                if key not in DEGs_SPs_dict:
                    DEGs_SPs_dict[key] = 0
                if key not in DAPs_SPs_dict:
                    DAPs_SPs_dict[key] = 0
            
            genes_dict = {'MONFIT':np.array(list(Preds_SPs_dict.values())), 'DAPs':np.array(list(DAPs_SPs_dict.values())),
                          'DEGs':np.array(list(DEGs_SPs_dict.values())), 'BG':np.array(list(BG_SPs_dict.values()))}
            colors = {'MONFIT':'lightgreen', 'DEGs':'cornflowerblue', 'DAPs':'gold', 'BG':'lightgrey'}
            
            maxx = np.max(BG_SP_dist)
            ind = np.arange(maxx) +1  
            SP_df = pd.DataFrame.from_dict(genes_dict)
            SP_df = SP_df.set_axis(ind)
            
            colors = {1:'lightgreen', 2:'cornflowerblue', 3:'gold', 4:'lightgrey'}
            cases = list(genes_dict.keys())
            
            # pvals = [pvalmwu_PD_BG, pvalmwu_DEG_BG, pvalmwu_DAP_BG, pvalmwu_PD_DEG, pvalmwu_PD_DAP]
            pairs=[("MONFIT", "BG"),("DEGs", "BG"), ('MONFIT','DEGs'),  ("DAPs", "BG"),('MONFIT','DAPs')]
            
            groups_c = ['MONFIT'] * len(TGs_SP_dist[0]) + ['DAPs'] * len(TGs_SP_dist[2]) + ['DEGs'] * len(TGs_SP_dist[1]) + ['BG'] * len(BG_SP_dist)       
            SP_c = ab = itertools.chain(TGs_SP_dist[0], TGs_SP_dist[2], TGs_SP_dist[1], BG_SP_dist)
                 
            df_annots = pd.DataFrame({'GeneGroup':groups_c, 'ShortPath':SP_c})
            

            SP_df_1dec = SP_df.apply(lambda x: round(x, 1))
            fig, ax = survey(dict(zip(SP_df_1dec.columns, SP_df_1dec.values.T)), SP_df_1dec.index, df_annots, cases, pairs)
            plt.savefig(f'{sd}/SP_PINK1_PPI_DEGsDAPsPreds_HSB.jpg', dpi = 350, format='jpg',  bbox_inches="tight")
            plt.show()
            plt.close()
            
            
            # are the first neighbors significant
            Statistics.write('FIRST NEIGHBOURS\n')
           
            Statistics.write(f'Is the number of the first neighbours of a {network} subgraph that a gene set forms with {target_gene}?\n')
            Statistics.write('Test: Hypergeometric test\n')

            pval, fold = firstneighSignificance('MONFIT', G, target_gene, MONFITpreds)
            Statistics.write('MONFIT\n')
            Statistics.write(f'fold = {fold}\t')            
            Statistics.write(f'pvalue = {pval}\n\n')
            pval, fold = firstneighSignificance('DEGs', G, target_gene, DEGs)
            Statistics.write('DEGs\n')
            Statistics.write(f'fold = {fold}\t')            
            Statistics.write(f'pvalue = {pval}\n\n')
            pval, fold = firstneighSignificance('DAPs', G, target_gene, DAPs)
            Statistics.write('DAPs\n')
            Statistics.write(f'fold = {fold}\t')            
            Statistics.write(f'pvalue = {pval}\n\n')

            Statistics.close()
   
            
