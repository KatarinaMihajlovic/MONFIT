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
import os, pickle
import pandas as pd
from statannotations.Annotator import Annotator
import itertools


def ShortestPath2Gene(G, gene, TargetGenes):
    G2gene = nx.single_source_dijkstra(G, gene)
    # print(G2gene)
    G2gene_lens = G2gene[0]
    G2gene_lens.pop(gene)

    TargetGenes_SP = {key:G2gene_lens[key] for key in G2gene_lens if key in TargetGenes}
    BG_SP = {key:G2gene_lens[key] for key in G2gene_lens if key not in TargetGenes}
    TargetGenes_SP_dist = np.array([TargetGenes_SP[key] for key in TargetGenes_SP])
    BG_SP_dist = np.array([BG_SP[key] for key in BG_SP])
    return TargetGenes_SP_dist, BG_SP_dist, TargetGenes_SP, BG_SP

def ShortestPath2Gene_MULTI(G, gene, TargetGenes):
    G2gene = nx.single_source_dijkstra(G, gene)
    # print(G2gene)
    G2gene_lens = G2gene[0]
    G2gene_lens.pop(gene)

    results_union = set().union(*TargetGenes)

    TargetGenes_SPs = []
    TargetGenes_SP_dists = []
    for l in TargetGenes:
        TargetGenes_SP = {key:G2gene_lens[key] for key in G2gene_lens if key in l}
        TargetGenes_SPs.append(TargetGenes_SP)
        TargetGenes_SP_dist = np.array([TargetGenes_SP[key] for key in TargetGenes_SP])
        TargetGenes_SP_dists.append(TargetGenes_SP_dist)

    BG_SP = {key:G2gene_lens[key] for key in G2gene_lens if key not in results_union}
    BG_SP_dist = np.array([BG_SP[key] for key in BG_SP])
    return TargetGenes_SP_dists, BG_SP_dist, TargetGenes_SPs, BG_SP



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
        
    print(observed_property)
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
    # ax.set_title(f'Shortest path length to PINK1 in {network} network', fontsize=34)
    ax.set_xlabel('Relative gene count (%)', fontsize=26)

    ax.invert_yaxis()
    ax.set_xticks([])
    # ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        rects = ax.barh(labels, widths, left=starts, height=0.6,
                        label=colname, color=color)

        # r, g, b, _ = color
        text_color = 'white' if color == 'cornflowerblue' else 'black'
        labels_ratio = [round(float(v), 1) if v > threshold else "" for v in rects.datavalues]  
        # print(labels_ratio)
        # labels_ratio = [round(float(x), 1) for x in labels_ratio if not isinstance(x, str)]
        # print(labels_ratio)

        ax.bar_label(rects, labels = labels_ratio, label_type='center', color=text_color, fontsize = 20, weight='bold', padding=5)
    ax.legend(ncols=len(category_names), bbox_to_anchor=(-0.03, 1), loc='lower left', title="Shortest path length", fontsize=24, title_fontsize=26)
    plt.yticks(fontsize=26)
    
    # threshold = 1
    # for c in ax.containers:
    #     # Filter the labels
    #     labels = [v if v > threshold else "" for v in c.datavalues]    
    #     ax.bar_label(rects, label_type='center', color=text_color, fontsize = 22, weight='bold', padding=10)

    
    annotator = Annotator(ax, pairs, data=df_annots, x='ShortPath', y='GeneGroup', order=order, orient='h')
    annotator.configure(test='Mann-Whitney-ls', text_format='star',fontsize=22, loc='outside')
    annotator.configure(comparisons_correction="BH", correction_format="replace")
    annotator.apply_and_annotate()
    return fig, ax

         
'''MAIN CODE'''

with open('input/CorePreds.pkl', 'rb') as handle:
    CorePreds = pickle.load(handle)
# with open('input/DEGs_Skupin.txt') as f:
#     DEGs_Skupin = f.read().splitlines()  
    

days_set = 'D8_D18_D25_D32_D37'
days = days_set.split('_')
days = [x[1:] for x in days]
    
All_genes = []
for file in os.listdir('input/Genelists'):
    day = file.split('_')[2][:-4]
    if day in days:
        Geneslist = pd.read_csv(f'input/Genelists/{file}')
        Geneslist = Geneslist['genes'].tolist()
        All_genes.append(Geneslist)
All_genes_inters = list(set.intersection(*map(set,All_genes)))
All_genes_Union = list(set.union(*map(set,All_genes)))



network = 'PPI'   
G = nx.read_edgelist('input/PPI_General_Biogrid_GeneSym.edgelist', data = False)
G  = G.subgraph(All_genes_Union)


target_gene = 'PINK1'

Statistics = open('output/Statistics.txt','w')
 
CorePreds_SP_dist, BG_SP_dist, CorePreds_SP, BG_SP = ShortestPath2Gene(G, target_gene, CorePreds)

# are the first neighbors significant
first_neigh = [key for key in CorePreds_SP if CorePreds_SP[key] == 1]
perc1st_neigh = len(first_neigh)/len(CorePreds)*100
write_txt(first_neigh, sd='output', save_file=f'{target_gene}_1neighCorePDPreds_{network}')
first_neigh_BG = [key for key in BG_SP if BG_SP[key] == 1]

X = len(first_neigh)
N = len(CorePreds)
K = len(first_neigh_BG) + len(first_neigh)
M = len(G.nodes())
fold = (X/N)/(K/M)
pval = hypergeom.sf(X-1, M, K, N)
print(f'1st Neighbours to {target_gene} in {network} significance: fold = {fold}, pval = {pval}')

Statistics.write(f'{network}\n') 
Statistics.write(f'1st Neighbors to {target_gene}, n = {len(first_neigh)} ({perc1st_neigh}%)\n')
Statistics.write(f'Are 1st Neighbors of {target_gene} from Core PD predictions significant?\n')
Statistics.write('Test: Hypergeometric test\n')
Statistics.write(f'fold = {fold}, pval = {pval}\t SIGNIFICANT\n')

# is the subnetwork more connected than expected at random 
prop = 'density'
p_value = SampWithReplac_TestSubgraph(G, CorePreds, target_gene, prop)
print(f'{prop} of Core PD preds with {target_gene} than expected at random ({network}): pvalue = {p_value}')
Statistics.write(f'Is the {prop} of Core PD preds with {target_gene} in {network} higher than expected at random?\n')
Statistics.write('Test: Sampling with replacement\n')
Statistics.write(f'pvalue = {p_value}\n\n')
  
subgraph_nodes_gene = CorePreds + [target_gene]
dens = nx.density(G.subgraph(subgraph_nodes_gene))


# MONFIT preds
MONFIT_preds = pd.read_csv('input/MO_GM_D8_D18_D25_D32_D37_1.533_PDpreds.csv', index_col=0)
MONFIT_preds = list(MONFIT_preds.index)
subgraph_nodes_gene = MONFIT_preds + [target_gene]
dens = nx.density(G.subgraph(subgraph_nodes_gene))



# ######### boxplot 


TGs_Skupin_SP_dist, BG_SP_dist, TGs_Skupin_SP, BG_SP = ShortestPath2Gene_MULTI(G, target_gene, [CorePreds, MONFIT_preds])

avgs = [np.mean(x) for x in TGs_Skupin_SP_dist]
print(avgs)

Preds_SP_cnt = np.bincount(TGs_Skupin_SP_dist[0])/len(TGs_Skupin_SP_dist[0])*100
BG_SP_cnt = np.bincount(BG_SP_dist)/len(BG_SP_dist)*100
MONFIT_preds_SP_cnt = np.bincount(TGs_Skupin_SP_dist[1])/len(TGs_Skupin_SP_dist[1])*100


Preds_SP_cnt = np.delete(Preds_SP_cnt, 0)
MONFIT_preds_SP_cnt = np.delete(MONFIT_preds_SP_cnt, 0)
BG_SP_cnt = np.delete(BG_SP_cnt, 0)

Preds_SPs = np.unique(TGs_Skupin_SP_dist[0])
MONFIT_preds_SPs = np.unique(TGs_Skupin_SP_dist[1])
BG_SPs = np.unique(BG_SP_dist)


### stacked bar chart

Preds_SPs_dict = dict(zip(Preds_SPs, Preds_SP_cnt))
MONFIT_preds_SPs_dict = dict(zip(MONFIT_preds_SPs, MONFIT_preds_SP_cnt))
BG_SPs_dict = dict(zip(BG_SPs, BG_SP_cnt))
for key in BG_SPs_dict:
    if key not in Preds_SPs_dict:
        Preds_SPs_dict[key] = 0
    if key not in MONFIT_preds_SPs_dict:
        MONFIT_preds_SPs_dict[key] = 0

genes_dict = {'MONFIT':np.array(list(MONFIT_preds_SPs_dict.values())), 'CorePDpreds':np.array(list(Preds_SPs_dict.values())),
             'BG':np.array(list(BG_SPs_dict.values()))}
colors = {'MONFIT':'lightgreen', 'CorePDpreds':'cornflowerblue', 'BG':'lightgrey'}

maxx = np.max(BG_SP_dist)
ind = np.arange(maxx) +1  
SP_df = pd.DataFrame.from_dict(genes_dict)
SP_df = SP_df.set_axis(ind)

colors = {1:'lightgreen', 2:'cornflowerblue', 3:'gold', 4:'lightgrey'}
cases = list(genes_dict.keys())

# pvals = [pvalmwu_PD_BG, pvalmwu_DEG_BG, pvalmwu_DAP_BG, pvalmwu_PD_DEG, pvalmwu_PD_DAP]
pairs=[("CorePDpreds", "BG"),("MONFIT", "BG"), ('MONFIT','CorePDpreds')]

groups_c = ['MONFIT'] * len(TGs_Skupin_SP_dist[1])+ ['CorePDpreds'] * len(TGs_Skupin_SP_dist[0])  + ['BG'] * len(BG_SP_dist)       
SP_c = itertools.chain(TGs_Skupin_SP_dist[1], TGs_Skupin_SP_dist[0], BG_SP_dist)
     
df_annots = pd.DataFrame({'GeneGroup':groups_c, 'ShortPath':SP_c})



SP_df_1dec = SP_df.apply(lambda x: round(x, 1))
fig, ax = survey(dict(zip(SP_df_1dec.columns, SP_df_1dec.values.T)), SP_df_1dec.index, df_annots, cases, pairs)
plt.savefig('output/SP_PINK1_PPI_MONFITPreds_HSB.jpg', dpi = 350, format='jpg',  bbox_inches="tight")
plt.show()
plt.close()

