# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:54:36 2024

@author: Katarina
"""

import os, pickle, math
from scipy.stats import hypergeom
import numpy as np
import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt

def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p
 

def load_GO(fname, genes):
    geneset=set(genes)
    GO_counter = {}
    gene2GO = {}
    #loading raw gene & go relationships
    fRead = open(fname, "r")
    for line in fRead.readlines():
        lspt = line.strip().split('\t')
        # print(lspt)
        if len(lspt)>1:
            try:
                geneid = lspt[0]
                term = lspt[1]
                if geneid in geneset:
                    if term not in GO_counter:
                        GO_counter[term] = 0
                    GO_counter[term] += 1
                    if geneid not in gene2GO:
                        gene2GO[geneid] = set()
                    gene2GO[geneid].add(term)
            except KeyError:
                pass
# 				print(lspt[0])
    fRead.close()
	#print(GO_counter)
	#filter out GO terms that annotates only one gene
    GO_counter2 = {}
    gene2GO2 = {}
    removed = set()
    for term in GO_counter:
        if GO_counter[term] > 1:
            GO_counter2[term] = GO_counter[term]
        else:
            removed.add(term)
    for gene in gene2GO:
        genego = gene2GO[gene].difference(removed)
        if len(genego)>0:
            gene2GO2[gene]=genego
    return [GO_counter2, gene2GO2]



def go_enrichment(clusters, gene2go, go_counts):
	M = len(gene2go) # total number of annotated genes in the dataset
	data = []
	i = -1
	enrichment = [[] for j in range(len(clusters))]
	cts = 0
	
	clts = []
	gos = []
	pvals = []
	
	NE_list = []
	
	for cluster in clusters:
		i+=1
		annotated_genes = []
		for gene in cluster:
			if gene in gene2go:
				annotated_genes.append(gene)
		N = len(annotated_genes) # number of annotated genes in the cluster
		#print N
		annotation_set = set()
		for gene in annotated_genes:
			for term in gene2go[gene]:
				annotation_set.add(term)
		#print len(annotation_set)
		for term in annotation_set:
			K = go_counts[term] # number of genes annotated with the given term in all the data
			X = 0	# number of gene in the clusters that are annotated with the go term
			for gene in annotated_genes:
				if term in gene2go[gene]:
					X += 1
			pval = hypergeom.sf(X-1, M, K, N) # faster and more accurate than 1-cdf
			clts.append(i)
			gos.append(term)
			pvals.append(pval)
			if pval <= 0.05:
				cts += 1
				#print "Cluster %i, term %s: X=%i, K=%i, pval = %s"%(i,term,X,K,str(pval))
				enrichment[i].append([term,pval])
			#print "%s %s"%(term,prb)

		#print(len(enrichment))    
		#Mean normalized entropy:
		d = float(len(annotation_set))	# number of different annotations in cluster c
		nec = 0.
		if d>1.:
			Cl = float(len(cluster))	# number of gene in cluster c
			sum_pi = 0.
			#counting map
			counting_map = {}
			for gene in cluster:
				if gene in gene2go:
					for go in gene2go[gene]:
						if go not in counting_map:
							counting_map[go] = 0.
						counting_map[go] += 1.
			for go in counting_map.keys():
				pi = counting_map[go]/Cl	# percentage of genes in c annotated with the considered go term
				sum_pi += pi*math.log(pi)
			nec = (-1./(math.log(d)))*sum_pi
		NE_list.append( nec )
	
	#applying BH correction (Benjamini-Hochner correction)
	BHpvals = p_adjust_bh(pvals)

	BHcts = 0
	BH_enrichment = [[] for j in range(len(clusters))]
	enriched_genes = []
	for i in range(len(clts)):
		#print pvals[i], BHpvals[i]
		if BHpvals[i]<=0.05:
			BHcts += 1
			BH_enrichment[clts[i]].append([gos[i],BHpvals[i]])
	for i in range(len(clusters)):
			cluster_set = set()
			enriched_gos = BH_enrichment[i]
			for gene in clusters[i]:
				for go, pval in enriched_gos:
					if gene in gene2go:
						if go in gene2go[gene]:
							cluster_set.add(gene)
			enriched_genes.append(list(cluster_set)) #genes that are enriched in each cluster
	
	#print cts, BHcts
	MNE = sum(NE_list)/float(len(NE_list))
	#print "Mean Normalized Entropy = %s"%(str(MNE))
	#print(BH_enrichment)
	enr_cluster = 0
	total_cluster = 0
	for i in range(len(clusters)):
		if len(clusters[i])>0:
			total_cluster += 1.
			if len(BH_enrichment[i])>0:
				enr_cluster += 1
	perc_enr_cluster = 100.*enr_cluster/total_cluster
	perc_genes = sum([len(enriched_genes[i]) for i in range(len(clusters))]) #number of genes enrihced in all clusters together (sum number of genes for each individual cluster)
	#print(perc_genes)
	perc_genes = 100.*perc_genes/float(len(gene2go))
	#print(perc_genes, float(len(gene2go)))    
	return [BH_enrichment, MNE, perc_genes, perc_enr_cluster]

def stats_sort_df(df, n, sort_by):
    # df['mean'] = df.mean(axis=1)
    cols = df.columns
    df = std_equivalent(df)
    df['mean'] = df[cols].mean(axis=1)
    df = df.sort_values(by=sort_by, ascending=False)
    df_n = df.head(n)
    return df, df_n    

def boxplot_df(df, title, score_type, save_file, tick_labels_sup=0):
    # boxplot
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.set_title(title, fontsize = 36)
    ax.set_ylabel(score_type, fontsize = 30)#, pad = 24)
    ax = sns.boxplot(data=df,  palette=sns.color_palette("husl", len(df.columns)))
    
    if tick_labels_sup == 0:
        tick_labels_sup = []
        for tick in df.columns:
            splt = tick.split('_')
            if len(splt) > 1:
                tick_sup = f'{splt[0]}$_{{{splt[1]}}}$'
                tick_labels_sup.append(tick_sup)
            else:
                tick_labels_sup.append(splt[0])
    
    ind = np.arange(len(tick_labels_sup))
    ax.set_xticks(ind)           
    ax.set_xticklabels(tick_labels_sup,  fontsize = 32, rotation=90)

    ax.tick_params(axis='x', which='major', labelsize=30)
    ax.tick_params(axis='y', which='major', labelsize=30)
    plt.subplots_adjust(top=1.1)
    plt.savefig(save_file, dpi = 350, format='jpg', bbox_inches='tight')	
    plt.show()
    plt.close() 

def std_equivalent(df): #not done
    medians = []
    min_errors = []
    max_errors = []
    for i in range(len(df.index)):
        median_val = np.percentile(df.iloc[i], 50)
        medians.append(median_val)
        q3, q1 = np.percentile(df.iloc[i], [84.13, 15.87])
        min_error = median_val - q1
        max_error = q3 - median_val 
        min_errors.append(min_error)
        max_errors.append(max_error)

    df['median'] = medians
    df['min_error'] = min_errors
    df['max_error'] = max_errors
    return df

 
def sns_heatmap(df_data, title, save_fig, hierarchical = False):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.tick_params(axis='x', which='major', labelsize=16)
    ax.tick_params(axis='y', which='major', labelsize=12)
    if hierarchical == False:       
        sns.heatmap(df_data, annot=True, cmap = 'coolwarm', fmt=".2f")#, vmin=0.3, vmax=0.7)
    else:
        sns.clustermap(df_data, annot=True, cmap = 'coolwarm', fmt=".2f")#, vmin=0.3, vmax=0.7)
    plt.title(title, fontsize = 24, pad = 24)
    plt.savefig(save_fig, dpi = 350, format='jpg', bbox_inches='tight')	
    plt.show()
    plt.close()

    
'''MAIN CODE'''

cell_conds = ['ND_8', 'ND_18', 'ND_25', 'ND_32', 'ND_37', 'WT_8', 'WT_18', 'WT_25', 'WT_32', 'WT_37']
PDpaths = 'PDmapPaths_noSBUK'


### ENRICHMENT IN DIFFERENT PD Pathways


sd = f'output/{PDpaths}'
try:
    os.makedirs(sd)
except FileExistsError:
    # directory already exists
    pass

for cc in cell_conds:
    try:
        genes = pd.read_csv(f'input/Genelists/Geneslist_{cc}.csv', header = 0)
        PPIlist = genes['genes'].to_list()
              
        fgoname = f'input/{PDpaths}.lst'
        go_counts, gene2go = load_GO(fgoname, PPIlist)
        print(f'{len(gene2go)} {PDpaths} go_counts genes')
        
        
        for wc in os.listdir(f'output/{cc}'):
            sd2 = f'{sd}/{cc}/{wc}'
            try:
                os.makedirs(sd2)
            except FileExistsError:
                # directory already exists
                pass
            
            if len(os.listdir(sd2)) != 10:
                j=0
                for file in os.listdir(f'output/{cc}/{wc}'):  
                    print(cc, file)
                    j+=1                
                    with open(f'output/{cc}/{wc}/{file}', 'rb') as handle:
                        G1_clust = pickle.load(handle)
                    Cp = [[str(y) for y in x] for x in G1_clust]           
            
                    bh, mne, eg, ec = go_enrichment(Cp, gene2go, go_counts)
                    terms = [[i[0] for i in nested] for nested in bh]
                    
                    
                    with open(f'{sd2}/{j}_ECEG.pkl', 'wb') as handle:
                        pickle.dump([ec,eg], handle)
    except Exception:
        pass


### PLOTTING

# PATHWAYS
n = 10       
sort_by = 'median'

PDpaths_l = ['PDmapPaths_noSBUK']
PDpaths_d = {'PDmapPaths_noSBUK':'$PD_{map}$'}
PDpaths_sub = {'EG_PDmapPaths_noSBUK':'$EG_{PD_{map}}$','EC_PDmapPaths_noSBUK':'$EC_{PD_{map}}$'}

cell_conds = ['ND_8', 'ND_18', 'ND_25', 'ND_32', 'ND_37', 'WT_8', 'WT_18', 'WT_25', 'WT_32', 'WT_37']

EC_EG = []

for PDpaths in PDpaths_l:
    print(PDpaths)
    PDpaths_ec = []
    PDpaths_eg = []
    
    for cc in cell_conds:
        print(cc)
        outdir = f'output/{PDpaths}/{cc}'
        cc_ec = []
        cc_eg = []
        wcs = os.listdir(outdir)
        print(len(wcs))
        for wc in os.listdir(outdir):
           wc_ec = []
           wc_eg = []
           for file in os.listdir(f'{outdir}/{wc}'):
               with open(f'{outdir}/{wc}/{file}', 'rb') as handle:
                   ECEG = pickle.load(handle) 
                   # print(ECEG)
                   wc_ec.append(ECEG[0])
                   wc_eg.append(ECEG[1])
           wc_ec_avg = np.mean(wc_ec)
           wc_eg_avg = np.mean(wc_eg)
           
           cc_ec.append(wc_ec_avg)
           cc_eg.append(wc_eg_avg)
       
        cc_ec = pd.DataFrame(cc_ec, columns=[cc], index=wcs)
        cc_eg = pd.DataFrame(cc_eg, columns=[cc], index=wcs)
       
        PDpaths_ec.append(cc_ec)
        PDpaths_eg.append(cc_eg)
    
    PDpaths_ec = pd.concat(PDpaths_ec, axis=1)
    PDpaths_eg = pd.concat(PDpaths_eg, axis=1)

    PDpaths_ec, PDpaths_ec_n = stats_sort_df(PDpaths_ec, n, sort_by)
    PDpaths_eg, PDpaths_eg_n = stats_sort_df(PDpaths_eg, n, sort_by)
    
    # EC.append(pd.DataFrame(PDpaths_ec[sort_by], columns=[PDpaths]))
    # EG.append(pd.DataFrame(PDpaths_eg[sort_by], columns=[PDpaths]))
    
    boxplot_df(PDpaths_ec, PDpaths_d[PDpaths], '% of enriched clusters', f'output/{PDpaths}/{PDpaths}_EC.jpg')
    boxplot_df(PDpaths_eg, PDpaths_d[PDpaths], '% of enriched genes', f'output/{PDpaths}/{PDpaths}_EG.jpg')

    sns_heatmap(PDpaths_ec_n, f'EC: {PDpaths_d[PDpaths]} - top {n} ({sort_by})', f'output/{PDpaths}/{PDpaths}_EC_top{n}_{sort_by}.jpg', hierarchical = False)    
    sns_heatmap(PDpaths_eg_n, f'EG: {PDpaths_d[PDpaths]} - top {n} ({sort_by})', f'output/{PDpaths}/{PDpaths}_EG_top{n}_{sort_by}.jpg', hierarchical = False)

    # metric = PDpaths_d[PDpaths]
    EC_EG.append(pd.DataFrame({f'EC_{PDpaths}': PDpaths_ec[sort_by]}))
    EC_EG.append(pd.DataFrame({f'EG_{PDpaths}': PDpaths_eg[sort_by]}))
   
    # EC_EG.append(pd.DataFrame(PDpaths_ec[sort_by], columns=[f'EC$_{{{metric}}}$']))
    # EC_EG.append(pd.DataFrame(PDpaths_eg[sort_by], columns=[f'EG$_{{{metric}}}$']))
    
    
EC_EG_df = pd.concat(EC_EG, axis=1)
# EC_EG_df.to_csv('output/Enrichment_WeightAVG.csv')

tick_labels = [PDpaths_sub[x] for x in EC_EG_df.columns]
boxplot_df(EC_EG_df, '', '% of enrichment', 'output/Enrichment_WeightAVG.jpg', tick_labels)

EC_EG_df.to_csv('output/Enrichment_WeightAVG.csv')

