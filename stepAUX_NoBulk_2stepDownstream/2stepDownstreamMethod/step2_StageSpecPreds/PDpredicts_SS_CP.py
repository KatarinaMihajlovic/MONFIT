# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 14:09:59 2022

@author: kmihajlo
"""
import os, pickle
import numpy as np
import pandas as pd
from scipy.spatial import distance
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kstest
from matplotlib.patches import Patch
from random import sample
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
 
with open('input/PDgenes_DGN_ALL.pkl', 'rb') as handle:
    DisGeNet_PDgenes = pickle.load(handle) 
    

def write_txt(list_l, file_path, name):
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    with open(f'{file_path}/{name}', 'w') as f: 
        for x in list_l:
            f.write(f'{x}\n')
            
def savePDPreds_txt(cc1cc2_genes, out_filename, save_dir = 'output/PD_Predictions'):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)        
    with open(f'{save_dir}/{out_filename}.txt', 'w') as outfile:
        for gene in cc1cc2_genes:
            outfile.write(f'{gene}\n')#f'{gene}\t{LitValid_AllGenes[gene][0]}\t{LitValid_AllGenes[gene][1]}\t{LitValid_AllGenes[gene][2]}\n')

def saveGenes_pkl(cc1cc2_genes, ccs2compare, save_dir, genes_type = 'All_genes'):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir) 
    with open(f'{save_dir}/{ccs2compare}_{genes_type}.pkl', 'wb') as fp:   #Pickling
        pickle.dump(cc1cc2_genes, fp)  

def gene_distances(CosSim_PD_vals, CosSim_C_vals, genes, day):
    dists = []
    for i in range(len(CosSim_PD_vals)):
        dist_PD = CosSim_PD_vals[i]
        dist_C = CosSim_C_vals[i]
        dist_gene = distance.euclidean(dist_PD, dist_C)
        dists.append(dist_gene)
    
    dists_pd = pd.DataFrame(dists, index = genes, columns = ['dist'])
    dists_pd = dists_pd.sort_values(by='dist', ascending=False)
 
    return dists_pd

def PDGs2Background_dist(dist_df, ccs2compare, cc1cc2_genes, stat_meas = 'MWU', PD_Gs = 'PD preds', DisGeNet_PDgenes = DisGeNet_PDgenes):
    print(PD_Gs)
    save_dir = 'output/StageSpecPreds'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)  
        
    Dist_gene2gene = []
    cc1cc2_genes_geneDist = []
    for gene in dist_df.index:
        if gene in cc1cc2_genes:
            cc1cc2_genes_geneDist.append(float(dist_df.loc[gene]))
        else:
            if PD_Gs == 'PD preds':
                if gene not in DisGeNet_PDgenes:
                    Dist_gene2gene.append(float(dist_df.loc[gene]))
            else:
                Dist_gene2gene.append(float(dist_df.loc[gene])) 
    
    # print(len(cc1cc2_genes_geneDist), len(Dist_gene2gene))
    if stat_meas == 'MWU':
        statmwu,pvalmwu = mannwhitneyu(cc1cc2_genes_geneDist, Dist_gene2gene, alternative='greater')
    elif stat_meas == 'KST':
        statmwu,pvalmwu = kstest(cc1cc2_genes_geneDist, Dist_gene2gene, alternative='greater')

    fold_dif = np.mean(cc1cc2_genes_geneDist)/np.mean(Dist_gene2gene)
    print(fold_dif)
    # print(pvalmwu)
    # if pvalmwu <= 0.05:
    #     print('YAAAY')
    ##computing the histograms
    num_bin = 50  
    lst = list(Dist_gene2gene)
    minv = min(lst)
    maxv = max(lst)
        
    bin_lims = np.linspace(minv,maxv,num_bin+1)
    bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    bin_widths = bin_lims[1:]-bin_lims[:-1]

    hist1, _ = np.histogram(cc1cc2_genes_geneDist, bins=bin_lims)
    hist2, _ = np.histogram(Dist_gene2gene, bins=bin_lims)
    print(np.sum(hist1), np.sum(hist2))
    ##normalizing
    hist1b = hist1/np.sum(hist1)*100 # divide with sum, the distribution surface sums to 1
    hist2b = hist2/np.sum(hist2)*100

    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_facecolor('xkcd:white')
    # ax.grid(color='grey', axis='x', linestyle='-', linewidth=0.25, alpha=0.5)
    # ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
    title = ccs2compare.split('D')
    day = title[1]
    # title = title[0] + '-P' + title[1]
    title = f'PD$_{{{day}}}$'

    ax.set_title(title, fontsize=40, pad=20) #': {PD_Gs} vs Other genes'   

    # sns.kdeplot(x=hist2b, multiple="stack", color='cornflowerblue', alpha=0.5)
    # sns.kdeplot(x=hist1b, multiple="stack", color='red', alpha = 0.5)

    ax.bar(bin_centers, hist1b, width = bin_widths, align = 'center', alpha = 0.5, color='red')
    ax.bar(bin_centers, hist2b, width = bin_widths, align = 'center', alpha = 0.5, color='cornflowerblue')

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    ax.set_xlabel('gene movement', fontsize=30, fontweight = 'bold') #dist_measure_l[dist_meas]
    ax.set_ylabel('relative count (%)', fontsize=30, fontweight = 'bold')
    for label in ax.get_xaxis().get_ticklabels()[::2]:
        label.set_visible(False)
    
    
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu <= 0.05:
        if stat_meas == 'MWU':
            fig.text(0.5, 0.52, f'MWU (1) > (2)\np-value = {pvalmwu_s}*\n', fontsize=28)
        elif stat_meas == 'KST':
            fig.text(0.5, 0.52, f'KST (1) > (2)\np-value = {pvalmwu_s}*\n', fontsize=28)
    else:
        if stat_meas == 'MWU':
            fig.text(0.5, 0.52, f'MWU (1) > (2)\np-value = {pvalmwu_s}\n', fontsize=28)
        elif stat_meas == 'KST':
            fig.text(0.5, 0.52, f'KST (1) > (2)\np-value = {pvalmwu_s}\n', fontsize=28)
    
    # if PD_Gs == 'PD preds':
    #     backg_genes = 'background'
    # else:
    #     backg_genes = 'Other genes'
    backg_genes = 'background'
    
    ax.legend(handles= [Patch(facecolor='red', alpha = 0.5, edgecolor='red',label=f'(1) {PD_Gs} ({len(cc1cc2_genes_geneDist)})'),
                Patch(facecolor='cornflowerblue', alpha = 0.5, edgecolor='cornflowerblue',label=f'(2) {backg_genes} ({len(Dist_gene2gene)})')],
              loc=1, labelspacing=1, prop={'size': 28}, facecolor ='white')
    
    plt.savefig(f'{save_dir}/{ccs2compare}_{PD_Gs}vsOG_{stat_meas}.jpg', dpi = 350, format='jpg')                        
    plt.show()
    plt.close()
    
def arccos_sim(cos_sim_all):
    arc_cos = np.arccos(cos_sim_all)
    # print(arc_cos)
    arc_cos = np.nan_to_num(arc_cos, copy=True, nan=0.0)
    arc_cos = 1 - arc_cos/np.pi
    return arc_cos
def squareroot_mf(cos_sim_all):
    sr = (1 - cos_sim_all)/2
    sr = np.sqrt(sr)
    # print(sr)
    sr = np.nan_to_num(sr, copy=True, nan=0.0)
    sr = 1 - sr
    return sr


in_dir = 'output/EuclideanDistance'


days = set()
for file in os.listdir(in_dir):
    day = file.split('_')[1]
    days.add(day)
days = list(days)


# Normalized Euclidean distance - not called taht, use instead eculdean distance normalized by dividing with norm

for day in days: 
    print(day)
    ccs2compare = f'WT{day}ND{day}'

    EucDist_PD = np.load(f'{in_dir}/ND_{day}_PWEucDist.npy')
    with open(f'output/ND_{day}_Genes.pkl', 'rb') as handle:
        Genes_PD = pickle.load(handle)  
    EucDist_PD_df = pd.DataFrame(EucDist_PD, index = Genes_PD, columns = Genes_PD)
    
    EucDist_C= np.load(f'{in_dir}/WT_{day}_PWEucDist.npy')
    with open(f'output/WT_{day}_Genes.pkl', 'rb') as handle:
        Genes_C = pickle.load(handle)  
    EucDist_C_df = pd.DataFrame(EucDist_C, index = Genes_C, columns = Genes_C)
   
    Genes_cc1cc2 = list(set(Genes_C) & set(Genes_PD))      
    saveGenes_pkl(Genes_cc1cc2, ccs2compare, save_dir = 'output/All_genes', genes_type = 'All_genes')   
    
    EucDist_PD_df = EucDist_PD_df[Genes_cc1cc2]
    EucDist_PD_df = EucDist_PD_df.loc[Genes_cc1cc2]
    EucDist_C_df = EucDist_C_df[Genes_cc1cc2]
    EucDist_C_df = EucDist_C_df.loc[Genes_cc1cc2]   
    
    sim_valsPD = EucDist_PD_vals = EucDist_PD_df.values
    sim_valsC = EucDist_C_vals = EucDist_C_df.values
    
    # Stage Specific Predictions

    with open(f'input/InitPreds/ND_{day}_InitPreds.pkl', 'rb') as handle:
        PD_preds = pickle.load(handle)  
               
    SS_PDpreds = list(set(Genes_C) & set(Genes_PD) & set(PD_preds))
    saveGenes_pkl(SS_PDpreds, ccs2compare, save_dir = 'output/StageSpecPreds', genes_type = 'PDpreds')
    savePDPreds_txt(SS_PDpreds, out_filename = f'{ccs2compare}_PDpreds', save_dir = 'output/StageSpecPreds')
    
    Removed_genes = list(set(PD_preds) - set(SS_PDpreds))
    write_txt(Removed_genes, file_path = 'output/StageSpecPreds/RemovedPDpreds', name = f'{ccs2compare}_RemPDpreds.txt')   

       
    dist_measure = gene_distances(sim_valsPD, sim_valsC, Genes_cc1cc2, day)
    
    for stat_meas in ['MWU']:#, 'KST']:
        ### check if PDpreds distances are significant compared to background
        PDGs2Background_dist(dist_measure, ccs2compare, SS_PDpreds, stat_meas = stat_meas, PD_Gs = 'PD predictions')
        ### check if PD known gene distances are significant compared to background
        PDGs2Background_dist(dist_measure, ccs2compare, DisGeNet_PDgenes, stat_meas = stat_meas, PD_Gs = 'DisGeNet PD genes')
    # save dist of SS_PDpreds  
    dist_measure = dist_measure.loc[SS_PDpreds].sort_values(['dist'], ascending=False)
    dist_measure.to_csv(f'output/StageSpecPreds/{ccs2compare}_PDpreds.csv', index = True, header = True)





 #### CORE PREDICTIONS    
    
PD_preds_ccs = []
PD_preds_ccs_dict = {}
for root, dirs, files in os.walk('output/StageSpecPreds'): 
    for file in files:
        if file.endswith('.pkl'):
            pair = file.split('_')[0]
            print(pair)
            with open(f'{root}/{file}', 'rb') as handle:
                cc_genes = pickle.load(handle)  
            PD_preds_ccs.append(cc_genes)
            PD_preds_ccs_dict[pair] = cc_genes
    CorePreds = list(set.intersection(*map(set,PD_preds_ccs)))
print(CorePreds)

   
# is the number of Core pd preds better than random  
reps = 10000
Successes = 0 
for rep in range(reps):
    PD_preds_ccs_rand = []  
    for pair in PD_preds_ccs_dict.keys():
        with open(f'output/All_genes/{pair}_All_genes.pkl', 'rb') as handle:
            cc1cc2_SGs = pickle.load(handle) 
        n_samp = len(PD_preds_ccs_dict[pair])
        lst = sample(cc1cc2_SGs, n_samp)
        PD_preds_ccs_rand.append(lst)
    CorePreds_rand = set.intersection(*map(set,PD_preds_ccs_rand))
    if len(CorePreds_rand) >= len(CorePreds):
        Successes+=1
pval = (Successes + 1)/(reps+1) 
print(pval)

# sort CorePreds by the avg distance across all stage pairs - from most distant to least
CorePreds_dist_pd = pd.DataFrame(index= CorePreds, columns = PD_preds_ccs_dict.keys())

for ccs2compare in PD_preds_ccs_dict.keys():
    print(ccs2compare)
    dist_df = pd.read_csv(f'output/StageSpecPreds/{ccs2compare}_PDpreds.csv', index_col = 0, header = 0)
    for gene in CorePreds:
        CorePreds_dist_pd.loc[gene,ccs2compare] = float(dist_df.loc[gene])
   
CorePreds_dist_pd['average'] = CorePreds_dist_pd.mean(axis=1)
CorePreds_dist_pd = CorePreds_dist_pd.sort_values(by='average', ascending=False)
CorePreds_dist_pd.to_csv(f'output/CorePreds.csv', index = True, header = True)



CorePreds_ordered = CorePreds_dist_pd.index.tolist()
with open('output/CorePreds.pkl', 'wb') as fp:   
    pickle.dump(CorePreds_ordered, fp)
savePDPreds_txt(CorePreds_ordered, out_filename = 'CorePreds', save_dir = 'output')



# comparing with MONFIT preds

PD_preds_GM = pd.read_csv('input/MO_GM_D8_D18_D25_D32_D37_1.157_PDpreds.csv', index_col=0)
PD_preds_MO = list(PD_preds_GM.index)

import venn

# Venn
genes = [set(CorePreds), set(PD_preds_MO)]
gene_sets = ['MONFIT', 'CorePreds']
labels = venn.generate_petal_labels(genes)#, fill=['number', 'logic'])
fig, ax = venn.venn2(labels, names=gene_sets, fontsize =24, figsize= (9, 12))
ax.set_title('Genes', fontsize = 28, y=0.9) # Gene Overlap between\nPD predictions, DEGs and DAPs'   
ax.get_legend().remove()
# leg = ax.legend(gene_sets, loc='lower left', bbox_to_anchor=(0, -0.14), fancybox=True, fontsize = 22)
# plt.savefig(f'{sd}/{save_file}_GeneOverlap_DEG{thresh_DEG}_DAP{thresh_DAPs}.jpg',  dpi = 350, format='jpg', bbox_inches='tight')	
plt.show()


print(list(set(PD_preds_MO)&set(CorePreds)))














# plt.figure()
# CorePreds_dist_pd.plot()