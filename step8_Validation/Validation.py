# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:12:57 2021

@author: kmihajlo
"""
import os, sys
import pickle
from scipy.stats import hypergeom, mannwhitneyu
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import pandas as pd
from statannotations.Annotator import Annotator

            
default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [99, 236, 243, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
]            
       

def Pubmedn(PDpreds_day_LitValid):
    cc1cc2_pubmed = []
    cc1cc2_genes_pubmed = []
    for gene in PDpreds_day_LitValid.keys():
        n = PDpreds_day_LitValid[gene][0]
        cc1cc2_pubmed.append(n)  
        cc1cc2_genes_pubmed.append([gene,n])
    return cc1cc2_pubmed, cc1cc2_genes_pubmed


def PDGs2Background_dist(Predictions_pubmed, BG_pubmed, title, save_dir):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)  
        
    statmwu,pvalmwu = mannwhitneyu(Predictions_pubmed, BG_pubmed, alternative='greater')

    ##computing the histograms
    num_bin = 50
    lst = list(Predictions_pubmed)
    minv = min(lst)
    maxv = max(lst)
        
    bin_lims = np.linspace(minv,maxv,num_bin+1)
    bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    bin_widths = bin_lims[1:]-bin_lims[:-1]

    hist1, _ = np.histogram(Predictions_pubmed, bins=bin_lims)
    hist2, _ = np.histogram(BG_pubmed, bins=bin_lims)
    print(np.sum(hist1), np.sum(hist2))
    ##normalizing
    hist1b = hist1/np.sum(hist1)*100 # divide with sum, the distribution surface sums to 1
    hist2b = hist2/np.sum(hist2)*100

    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_facecolor('xkcd:white')        
    # ax.set_title('Co-occurence of genes with\nParkinson\'s disease term in Pubmed', fontsize=34, pad=20) #': {PD_Gs} vs Other genes'   

    ax.bar(bin_centers, hist1b, width = bin_widths, align = 'center', alpha = 0.5, color='red')
    ax.bar(bin_centers, hist2b, width = bin_widths, align = 'center', alpha = 0.5, color='cornflowerblue')
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    ax.xaxis.offsetText.set_fontsize(20)
     
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax.set_xlabel('log(PubMed co-occurences + 1)', fontsize=26) #dist_measure_l[dist_meas]
    ax.set_ylabel('Relative count (%)', fontsize=26)
    for label in ax.get_xaxis().get_ticklabels()[::2]:
        label.set_visible(False)
    
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu <= 0.05:
        fig.text(0.5, 0.52, f'MWU (1) > (2)\np-value = {pvalmwu_s}*\n', fontsize=28)
    else:
        fig.text(0.5, 0.52, f'MWU (1) > (2)\np-value = {pvalmwu_s}\n', fontsize=28)
    
    ax.legend(handles= [Patch(facecolor='red', alpha = 0.5, edgecolor='red',label=f'(1) MONFIT predictions ({len(Predictions_pubmed)})'),
                Patch(facecolor='cornflowerblue', alpha = 0.5, edgecolor='cornflowerblue',label=f'(2) background ({len(BG_pubmed)})')],
              loc=1, labelspacing=1, prop={'size': 28}, facecolor ='white')
    
    plt.savefig(f'{save_dir}/{title}_PDpredsvsOG_MWU.jpg', dpi = 350, format='jpg',  bbox_inches="tight")        
                
    plt.show()
    plt.close()
 
    
    
def nvalid_genes(genes_LitValid):
    n_valids = 0
    for gene in genes_LitValid:
        if genes_LitValid[gene][0] > 0:
            n_valids +=1
        elif genes_LitValid[gene][1] != '':
            n_valids +=1
    return n_valids



def PDenrichedClust(gene_list, total_genes, PD_genes):
    Possible_PDgenes = list(set(total_genes) & set(PD_genes))                          
    PDgenes_clust = [x for x in Possible_PDgenes if x in gene_list]                            
    
    M = len(total_genes)
    K = len(Possible_PDgenes)
    N = len(gene_list)
    X = len(PDgenes_clust)
    # print(M, K, N, X)
    try:
        fold = (X/N)/(K/M)
    except ZeroDivisionError:
        fold = None
        percPDgenes = None
        pval = None
    if fold >= 1:
        # print(fold)
        pval = hypergeom.sf(X-1, M, K, N)
        if pval <= 0.05: 
            percPDgenes = len(PDgenes_clust)/len(Possible_PDgenes)*100
        else:
            percPDgenes = None
    else:
        pval = None
        percPDgenes = None
    return fold, percPDgenes, pval




def plotPubMed(cc1cc2_pubmed, cc1cc2_AllGenes_pubmed, cc_pair, sd):
    cc1cc2_pubmed = np.array(cc1cc2_pubmed)
    cc1cc2_AllGenes_pubmed = np.array(cc1cc2_AllGenes_pubmed)
    statmwu,pvalmwu = mannwhitneyu(cc1cc2_pubmed,cc1cc2_AllGenes_pubmed, alternative='greater')
    # print(pvalmwu)
      
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_title('Co-occurence of genes with Parkinson\'s disease term in Pubmed', fontsize=34, pad=20)    

    plt.hist(cc1cc2_pubmed, bins = 100, label=f'(1) PD predictions ({len(cc1cc2_pubmed)})', alpha = 0.5, color='red')#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 
    plt.hist(cc1cc2_AllGenes_pubmed, bins = 100, label=f'(2) background ({len(cc1cc2_AllGenes_pubmed)})', alpha = 0.5, color='cornflowerblue')#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 

    # plt.axvline(x = v_line, color = 'r', label = f'{TP}% thresh ({TP_len})')
    ax.set_yscale('log')
    ax.set_xscale('log')
   
    # plt.title(cc_pair, fontsize=26) 
    ax.set_xlabel('PubMed co-occurence', fontsize=26) #dist_measure_l[dist_meas]
    ax.set_ylabel('frequency', fontsize=26)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=26)
        
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu < 0.05:
        fig.text(0.5, 0.56, f'MWU (1) > (2)\np-value = {pvalmwu_s}*', fontsize=28)
    else:
        fig.text(0.5, 0.56, f'MWU (1) > (2)\np-value = {pvalmwu_s}', fontsize=28)

    plt.savefig(f'{sd}/{cc_pair}_PDpredsvsOGs.jpg', dpi = 350, format='jpg')                        
    plt.show()
    plt.close()
  
 

# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p


    
#plot precision and recall at different tops as different plots

####### MAIN CODE
wd = str(sys.argv[1])

with open(f'{wd}/input/PDgenes_PDmap.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)
with open(f'{wd}/input/Gene4PD.pkl', 'rb') as handle:
    Gene4PD = pickle.load(handle)
with open(f'{wd}/input/PDgenes_DGN_ALL.pkl', 'rb') as handle:
    DGN_PDgenes_ALL = pickle.load(handle)


PDgenes_union = list(set().union(Gene4PD, DGN_PDgenes_ALL))
    

with open(f'{wd}/input/LitValid_AllGenes.pkl', 'rb') as handle:
    LitValid_AllGenes = pickle.load(handle) 

# for gene in LitValid_AllGenes:
#     LitValid_AllGenes[gene][1] = ''
#     if gene in PD_genes:
#         LitValid_AllGenes[gene][1] = 'PD_map'
 
    
   
in_dir = f'{wd}/input/Predictions'

for days_set in os.listdir(in_dir):   
    # days_set = 'D18_D25_D32_D37'
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
    
    All_CommonGenes_LitValid = {}
    for key in All_genes_inters:
        try:
            All_CommonGenes_LitValid[key] = LitValid_AllGenes[key]
        except Exception:
            print(key)


    PDgenes_union_inters = list(set(PDgenes_union) & set(All_genes_inters))
    DGN_PDgenes_ALL_inters = list(set(DGN_PDgenes_ALL) & set(All_genes_inters))
    Gene4PD_inters = list(set(Gene4PD) & set(All_genes_inters))
    print(len(Gene4PD_inters), len(DGN_PDgenes_ALL_inters), len(PDgenes_union_inters))
   
    
    for perc_case in os.listdir(f'{in_dir}/{days_set}'):
        for file in os.listdir(f'{in_dir}/{days_set}/{perc_case}'):
            # perc_case = '1.476perc'
            # file = 'GM_Length_D18_D25_D32_D37_1.476_PDpreds.csv'
            save_file = file[:-4]          
            print(save_file)
            pvals = []
            PDgenes_cases = []
           
            sd = f'{wd}/output/{days_set}/{perc_case}'
            if not os.path.exists(sd):
                os.makedirs(sd) 
                
            # with open(f'{sd}/All_genes_inters.pkl', 'wb') as fp:   
            #     pickle.dump(All_genes_inters, fp)  

                
            S2Ss_df = pd.read_csv(f'{in_dir}/{days_set}/{perc_case}/{file}', header = 0, index_col = 0)
            S2Ss = S2Ss_df.index.tolist()
            S2Ss_LitValid = dict((k, LitValid_AllGenes[k]) for k in S2Ss if k in LitValid_AllGenes)                              
            PDgenes_size = []
            
            # % and enrichment of pd genes from Gene4PD
            Gene4PD_genes = [gene for gene in S2Ss if gene in Gene4PD_inters]
            PDgenes_perc1 = len(Gene4PD_genes)/len(S2Ss)*100
            fold1, percPDgenes1, pval1 = PDenrichedClust(S2Ss, All_genes_inters, Gene4PD_inters) #percPDgenes compared to all the PDgenes from a certain set of pd genes
            # PDgenes_enr_file.write(f'Gene4PD\n% of genes that are Gene4PD genes: {PDgenes_perc}\nfold = {fold}\npvalue = {pval}\n\n')
            pvals.append(pval1)
            PDgenes_cases.append('Gene4PD')
            PDgenes_size.append(len(Gene4PD_genes))
                    
            # % and enrichment of pd genes from Disgenet_ALL
            DGN_genes_ALL = [gene for gene in S2Ss if gene in DGN_PDgenes_ALL_inters]
            PDgenes_perc2 = len(DGN_genes_ALL)/len(S2Ss)*100
            fold2, percPDgenes2, pval2 = PDenrichedClust(S2Ss, All_genes_inters, DGN_PDgenes_ALL_inters) #percPDgenes compared to all the PDgenes from a certain set of pd genes
            # PDgenes_enr_file.write(f'DisGeNet_ALL\n% of genes that are DisGeNet genes: {PDgenes_perc}\nfold = {fold}\npvalue = {pval}\n\n')
            pvals.append(pval2)
            PDgenes_cases.append('DisGeNet')
            PDgenes_size.append(len(DGN_genes_ALL))


            # % of union of all the databases
            PDgenes_union_S2S = [gene for gene in S2Ss if gene in PDgenes_union_inters]
            print('PDgenes_union')
            PDgenes_perc3 = len(PDgenes_union_S2S)/len(S2Ss)*100
            print(PDgenes_perc3)
            fold3, percPDgenes3, pval3 = PDenrichedClust(S2Ss, All_genes_inters, PDgenes_union_inters) #percPDgenes compared to all the PDgenes from a certain set of pd genes
            # PDgenes_enr_file.write(f'Union of DBs\n% of genes that are PD genes: {PDgenes_perc}\nfold = {fold}\npvalue = {pval}\n\n')
            pvals.append(pval3)
            PDgenes_cases.append('Union')
            PDgenes_size.append(len(PDgenes_union_S2S))

            # adjust for mutiple-hypothesis testing
            padj = p_adjust_bh(pvals)

            PDgenes_enr_file = open(f'{sd}/{save_file}_PDgenesEnrich.txt','w')        
            PDgenes_enr_file.write(f'Gene4PD\n% of genes that are Gene4PD genes: {PDgenes_perc1}\nfold = {fold1}\npvalue = {padj[0]}\n\n')
            PDgenes_enr_file.write(f'DisGeNet_ALL\n% of genes that are DisGeNet genes: {PDgenes_perc2}\nfold = {fold2}\npvalue = {padj[1]}\n\n')
            PDgenes_enr_file.write(f'Union of DBs\n% of genes that are PD genes: {PDgenes_perc3}\nfold = {fold3}\npvalue = {padj[2]}\n\n')
            PDgenes_enr_file.close()
                
            
         

            
            # validation of genes
            BG_CorePreds = [gene for gene in All_genes_inters if gene not in S2Ss]
            # BG_CorePreds_LitValid = {your_key: LitValid_AllGenes[your_key] for your_key in BG_CorePreds}
            BG_CorePreds_LitValid = {}
            for key in BG_CorePreds:
                try:
                    BG_CorePreds_LitValid[key] = LitValid_AllGenes[key]
                except Exception:
                    print(key)


                
            with open(f'{sd}/{save_file}_LitValid.txt', 'w') as f:
                for gene in S2Ss:
                    try:
                        f.write(f'{gene}\t{S2Ss_LitValid[gene][0]}\n')
                    except:
                        f.write(f'{gene}\t nan \n')
            with open(f'{sd}/{save_file}_LitValid.pkl', 'wb') as fp:   
                pickle.dump(S2Ss_LitValid, fp)               
            with open(f'{sd}/All_CommonGenes_LitValid.pkl', 'wb') as fp:   
                pickle.dump(All_CommonGenes_LitValid, fp)
            
                 
            # simple hist
            S2Ss_pubmed,_ = Pubmedn(S2Ss_LitValid)
            BG_pubmed,_ = Pubmedn(BG_CorePreds_LitValid) #all genes a OtherGenes_LitValid
            gene_group = save_file[10:]
            # plotPubMed(S2Ss_pubmed, BG_pubmed, gene_group, sd)
            PDGs2Background_dist(S2Ss_pubmed, BG_pubmed, gene_group, sd)
            
            
            #  enrichment in known PD genes
            PDgenes_size.insert(0, len(S2Ss))
            PDgenes_cases.insert(0, 'MONFIT')
            pairs=[("MONFIT", "Gene4PD"), ("MONFIT", "DisGeNet"), ("MONFIT", "Union")]

            
          
            default_colors[5][3] = 0.45
            default_colors[4][3] = 0.45
            colors_c = [default_colors[0], default_colors[2], default_colors[5],default_colors[4]]

            df_Enr =  pd.DataFrame({'GeneGroup':PDgenes_cases, 'SizeGenes':PDgenes_size})
            fontsize_ls = 32
            
            fig, ax = plt.subplots(figsize=(8, 9))
            # Remove top and right border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_facecolor('xkcd:white')
            
            
            ax = df_Enr['SizeGenes'].plot.bar(x='GeneGroup', y='SizeGenes',  color=colors_c,  rot=0, width=1)
            ax.set_xticklabels(df_Enr['GeneGroup'], fontsize = 26, rotation=90)
            ax.set_ylabel('Number of genes', fontsize = fontsize_ls)
            
            ax.tick_params(axis='y', which='major', labelsize=26)
            plt.ylim(0,len(S2Ss)+10)
            
            annotator = Annotator(ax, pairs, data=df_Enr, x='GeneGroup', y='SizeGenes', order=PDgenes_cases,fontsize=32, loc='outside')
            
            annotator.configure(fontsize=32)

            new_annoations = []
            print(padj)
            for p in padj:
                if 5.00e-02 < p <= 1.00e+00:
                    new_annoations.append('ns')
                elif 1.00e-02 < p <= 5.00e-02: 
                    new_annoations.append('*')
                elif 1.00e-03 < p <= 1.00e-02: 
                    new_annoations.append('**')
                elif 1.00e-04 < p <= 1.00e-03: 
                    new_annoations.append('***')
                elif p <= 1.00e-04: 
                    new_annoations.append('****')
            annotator.annotate_custom_annotations(new_annoations)
            
            plt.savefig(f'{sd}/LitValid_Enrich.jpg',  dpi = 350, format='jpg', bbox_inches='tight')	
            plt.show()
            plt.close()
                     
         
        