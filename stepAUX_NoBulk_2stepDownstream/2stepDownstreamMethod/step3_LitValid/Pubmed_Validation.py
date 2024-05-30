# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:12:57 2021

@author: kmihajlo
"""
import os, math
import pickle
from scipy.stats import hypergeom, mannwhitneyu
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

def Genesym2Ensembl(filename = 'input/Homo_sapiens.gene_info'):
    Genesym2Ensembl = {}
    with open(filename) as f:
        f.readline()
        for line in f:
            lspt = line.split('\t')
            GeneNames = lspt[5].split('|')
            for Genename in GeneNames:
                if 'Ensembl' in Genename:
                    #print(Genename)
                    Ensembl = Genename.split(':')[1]
                    Genesym2Ensembl[lspt[2]] = Ensembl
    return Genesym2Ensembl

Genesym2Ensembl = Genesym2Ensembl()

def logPubmedn(PDpreds_day_LitValid):
    cc1cc2_pubmed = []
    cc1cc2_genes_pubmed = []
    for gene in PDpreds_day_LitValid.keys():
        logn = math.log(PDpreds_day_LitValid[gene][0] + 1, 10)
        cc1cc2_pubmed.append(logn)  
        cc1cc2_genes_pubmed.append([gene,logn])
    return cc1cc2_pubmed, cc1cc2_genes_pubmed


def Pubmedn(PDpreds_day_LitValid):
    cc1cc2_pubmed = []
    cc1cc2_genes_pubmed = []
    for gene in PDpreds_day_LitValid.keys():
        n = PDpreds_day_LitValid[gene][0]
        cc1cc2_pubmed.append(n)  
        cc1cc2_genes_pubmed.append([gene,n])
    return cc1cc2_pubmed, cc1cc2_genes_pubmed

 
def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)


    
def nvalid_genes(genes_LitValid):
    n_valids = 0
    for gene in genes_LitValid:
        if genes_LitValid[gene][0] > 0:
            n_valids +=1
        elif genes_LitValid[gene][1] != '':
            n_valids +=1
    return n_valids


def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)

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
        pval = 1
    if fold >= 1:
        # print(fold)
        pval = hypergeom.sf(X-1, M, K, N)
        if pval <= 0.05: 
            percPDgenes = len(PDgenes_clust)/len(Possible_PDgenes)*100
        else:
            percPDgenes = 0
    else:
        pval = 1
        percPDgenes = 0
    return fold, percPDgenes, pval


def LitEnrich(All_litValid, gene_list, outf, s2s):  # not working properly
    M = 0
    K = 0
    N = 0
    X = 0
    LitEnr_genes = []
    for gene in All_litValid:
        M+=1
        if gene in gene_list:
            N+=1
        if All_litValid[gene][0] > 0:
        # if All_litValid[gene][1] != '':
        # if All_litValid[gene][0] > 0:
            K+=1
            if gene in gene_list:
                X+=1
    
    print(M,K,N,X)
    outf.write(f'{M}\t{K}\t{N}\t{X}\n')
    
    perc = X/N*100
    print(perc)
    outf.write(f'{perc}%\n')
    
    try:
        fold = (X/N)/(K/M)
        print(fold)
        outf.write(f'{fold}\n')    
    except ZeroDivisionError:
        fold = 0
        pval = 1
    
    if fold >= 1:
        pval = hypergeom.sf(X-1, M, K, N)
        print(pval)
        outf.write(f'{pval}\n')
    
        if pval <= 0.05: 
            #print(f'Enriched clust {i},  fold = {fold},  pval = {pval}')
            #print(N, X)
            print('Enriched in LitValid - PubMed,PD_map')    
            outf.write('Enriched in LitValid - PubMed,PD_map\n')   
            print(f'percentage of enriched genes: {X/N*100}')
            print(f'percentage of enriched background: {K/M*100}')
            outf.write(f'percentage of enriched genes: {X/N*100}\n')
            outf.write(f'percentage of enriched background: {K/M*100}\n')
    else:
        pval = 1
            
    print('\n')  
    LitEnr_genes = [s2s, N, X, pval]
    return LitEnr_genes

colorp = {"Other_gene":'cornflowerblue', "PDpred_valid":'red'}
alphap = {"Other_gene":0.2, "PDpred_valid":0.8}
markerp = {"Other_gene":'o', "PDpred_valid":'x'}





def plotPubMed(cc1cc2_pubmed, cc1cc2_AllGenes_pubmed, cc_pair, sd):
    cc1cc2_pubmed = np.array(cc1cc2_pubmed)
    cc1cc2_AllGenes_pubmed = np.array(cc1cc2_AllGenes_pubmed)
    statmwu,pvalmwu = mannwhitneyu(cc1cc2_pubmed,cc1cc2_AllGenes_pubmed, alternative='greater')
    # print(pvalmwu)
      
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_title(cc_pair, fontsize=34, pad=20)    

    plt.hist(cc1cc2_pubmed, bins = 100, label=f'(1) PD predictions ({len(cc1cc2_pubmed)})', alpha = 0.5, color='red')#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 
    plt.hist(cc1cc2_AllGenes_pubmed, bins = 100, label=f'(2) background ({len(cc1cc2_AllGenes_pubmed)})', alpha = 0.5, color='cornflowerblue')#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 

    # plt.axvline(x = v_line, color = 'r', label = f'{TP}% thresh ({TP_len})')
    ax.set_yscale('log')
    # ax.set_xscale('log')
  
    plt.title(cc_pair, fontsize=26) 
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
    ax.set_title('Co-occurence of genes with\nParkinson\'s disease term in Pubmed', fontsize=34, pad=20) #': {PD_Gs} vs Other genes'   

    ax.bar(bin_centers, hist1b, width = bin_widths, align = 'center', alpha = 0.5, color='red')
    ax.bar(bin_centers, hist2b, width = bin_widths, align = 'center', alpha = 0.5, color='cornflowerblue')
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    ax.xaxis.offsetText.set_fontsize(20)
     
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax.set_xlabel('log(PubMed co-occurences + 1)', fontsize=26, fontweight = 'bold') #dist_measure_l[dist_meas]
    ax.set_ylabel('relative count (%)', fontsize=26, fontweight = 'bold')
    for label in ax.get_xaxis().get_ticklabels()[::2]:
        label.set_visible(False)
    
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu <= 0.05:
        fig.text(0.5, 0.52, f'MWU (1) > (2)\np-value = {pvalmwu_s}*\n', fontsize=28)
    else:
        fig.text(0.5, 0.52, f'MWU (1) > (2)\np-value = {pvalmwu_s}\n', fontsize=28)
    
    ax.legend(handles= [Patch(facecolor='red', alpha = 0.5, edgecolor='red',label=f'(1) PD predictions ({len(Predictions_pubmed)})'),
                Patch(facecolor='cornflowerblue', alpha = 0.5, edgecolor='cornflowerblue',label=f'(2) background ({len(BG_pubmed)})')],
              loc=1, labelspacing=1, prop={'size': 28}, facecolor ='white')
    
    plt.savefig(f'{save_dir}/{title}_PDpredsvsOG_MWU.jpg', dpi = 350, format='jpg',  bbox_inches="tight")        
                
    plt.show()
    plt.close()
 

    

####### MAIN CODE

with open('input/Gene4PD.pkl', 'rb') as handle:
    Gene4PD = pickle.load(handle)


with open('input/LitValid_AllGenes.pkl', 'rb') as handle:
    LitValid_AllGenes = pickle.load(handle) 


All_genes_inters = []
for file in os.listdir('input/Geneslist'):
    with open(f'input/Geneslist/{file}', 'rb') as handle:
        Geneslist = pickle.load(handle)        
    All_genes_inters.append(Geneslist)
All_genes_inters = list(set.intersection(*map(set,All_genes_inters)))


Gene4PD_inters = list(set(Gene4PD) & set(All_genes_inters))
   

with open('input/CorePreds.pkl', 'rb') as handle:
    S2Ss = pickle.load(handle)

S2Ss_LitValid = dict((k, LitValid_AllGenes[k]) for k in S2Ss if k in LitValid_AllGenes)                              
sd = 'output'
save_file = 'CorePreds'

PDgenes_size = []
pvals = []
PDgenes_cases = []
            
PDgenes_enr_file = open(f'{sd}/{save_file}_PDgenesEnrich.txt','w')
# % and enrichment of pd genes from Gene4PD
Gene4PD_genes = [gene for gene in S2Ss if gene in Gene4PD_inters]
PDgenes_perc = len(Gene4PD_genes)/len(S2Ss)*100
print('Gene4PD')
print(PDgenes_perc)
fold, percPDgenes, pval = PDenrichedClust(S2Ss, All_genes_inters, Gene4PD_inters) #percPDgenes compared to all the PDgenes from a certain set of pd genes
print(fold, pval, percPDgenes)
PDgenes_enr_file.write(f'Gene4PD\n% of genes that are Gene4PD genes: {PDgenes_perc}\nfold = {fold}\npvalue = {pval}\n\n')
pvals.append(pval)
PDgenes_cases.append('Gene4PD')
PDgenes_size.append(len(Gene4PD_genes))
        

PDgenes_enr_file.close()

S2S_list = [len(S2Ss)]*len(pvals)
# validation of genes
All_CommonGenes_LitValid = {your_key: LitValid_AllGenes[your_key] for your_key in All_genes_inters}
BG_CorePreds = [gene for gene in All_genes_inters if gene not in S2Ss]
BG_CorePreds_LitValid = {your_key: LitValid_AllGenes[your_key] for your_key in BG_CorePreds}


    
with open(f'{sd}/{save_file}_LitValid.txt', 'w') as f:
    for gene in S2Ss:
        f.write(f'{gene}\t{S2Ss_LitValid[gene][0]}\n')
with open(f'{sd}/{save_file}_LitValid.pkl', 'wb') as fp:   
    pickle.dump(S2Ss_LitValid, fp)               
with open(f'{sd}/All_CommonGenes_LitValid.pkl', 'wb') as fp:   
    pickle.dump(All_CommonGenes_LitValid, fp)

     

# simple hist
S2Ss_pubmed,_ = Pubmedn(S2Ss_LitValid)
BG_pubmed,_ = Pubmedn(BG_CorePreds_LitValid) #all genes a OtherGenes_LitValid
gene_group = save_file
# plotPubMed(S2Ss_pubmed, BG_pubmed, gene_group, sd)
PDGs2Background_dist(S2Ss_pubmed, BG_pubmed, gene_group, sd)

### check enrichment of PD preds with PD proof from Pubmed and PD_map, 
file_out = f'{save_file}_LitValid_Enirch.txt'
outf = open(f'{sd}/{file_out}','w')
LitEnr_genes = LitEnrich(All_CommonGenes_LitValid, S2Ss_LitValid, outf, save_file)
outf.close()        



# enrichment in known PD genes
width = 0.3    
opacity = 0.7
fontsize_ls = 32

fig, ax = plt.subplots(figsize=(16, 9))
# Remove top and right border
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
# Remove y-axis tick marks
ax.yaxis.set_ticks_position('none')
# Set plot title

ax.set_title('Enrichment of PD gene predictions in known PD genes', fontsize = 36)
# Add major gridlines in the y-axis
ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)
ax.set_facecolor('xkcd:white')

bar1 = ax.bar(np.arange(len(S2S_list)), S2S_list, width, align='center', alpha=opacity, color='r', label='PD gene predictions')
bar2 = ax.bar(np.arange(len(PDgenes_size)) + width, PDgenes_size, width, align='center', alpha=opacity,  color='limegreen', label='Known PD genes')
 
ax.legend(fontsize=30, loc = 'upper left')
ax.set_ylabel('#genes', fontsize = fontsize_ls, fontweight = 'bold')
ax.set_xticks(np.arange(len(S2S_list))+width/2)
ax.set_xticklabels(PDgenes_cases,  fontsize = fontsize_ls) 
ax.tick_params(axis='y', which='major', labelsize=fontsize_ls)
plt.ylim(0,550)
i = 0
for rect in bar1:
    height = rect.get_height()
    if float(pvals[i]) <= 0.05:
        pval = "{:.3e}".format(pvals[i])
        plt.text(rect.get_x() + rect.get_width(), height+0.1, f'{pval}*', ha='center', va='bottom', fontsize = fontsize_ls)
    else:
        plt.text(rect.get_x() + rect.get_width(), height+0.1, 'ns', ha='center', va='bottom', fontsize = fontsize_ls)
    i+=1

plt.savefig(f'{sd}/LitValid_Enrich.jpg',  dpi = 350, format='jpg', bbox_inches='tight')	
plt.show()
plt.close()
         
         
            
            
            
            










