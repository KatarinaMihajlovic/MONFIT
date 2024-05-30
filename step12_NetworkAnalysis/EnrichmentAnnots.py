# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:15:35 2023

@author: kimhajlo
"""

#import networkx as nx
import math
from scipy.stats import hypergeom
import numpy as np
import matplotlib.pyplot as plt
import os.path
import pandas as pd
from matplotlib_venn import venn2
import pickle

from goatools.obo_parser import GODag
godag = GODag(obo_file="input/go-basic.obo") 

with open('input/Entrez2Sym_name_synon.pkl', 'rb') as handle:
    Entrez2Sym = pickle.load(handle)
     

# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p

# Read gene ontology
def load_gene_ontology(fname):
	go2name = {}
	ifile = open(fname,'r')
	go = ""
	name = ""
	for line in ifile.readlines():
		try:
			if len(line)>3:
				if line[:3] == "id:":
					lspt = line.strip().split(':')
					go = "GO:"+lspt[2]
			if len(line)>5:
				if line[:5] == "name:":
					lspt = line.strip().split(':')
					name = lspt[1]
			if go != "" and name != "":
				go2name[go]=name
				go = ""
				name = ""
		except:
			print( line)
	ifile.close()
	#print "%i annotation loaded"%(len(go2name))
	return go2name


def load_PD(fname, genes):
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

# Read GO terms annotations
def load_GO(fname, genes):
	geneset=set(genes)
	GO_counter = {}
	gene2GO = {}
	#loading raw gene & go relationships
	fRead = open(fname, "r")
	for line in fRead.readlines():
		lspt = line.strip().split()
		if len(lspt)>1:
			geneid = lspt[0]
			term = lspt[1]
			try:
				geneid = Entrez2Sym[geneid]
			except KeyError:
				geneid = ''
			if geneid in geneset:
				if term not in GO_counter:
					GO_counter[term] = 0
				GO_counter[term] += 1
				if geneid not in gene2GO:
					gene2GO[geneid] = set()
				gene2GO[geneid].add(term)
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
	
	cls = []
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
			cls.append(i)
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
	for i in range(len(cls)):
		#print pvals[i], BHpvals[i]
		if BHpvals[i]<=0.05:
			BHcts += 1
			BH_enrichment[cls[i]].append([gos[i],BHpvals[i]])
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
# 	print(len(enriched_genes[0]), len(clusters))
# 	print(perc_genes)
    #print(perc_genes)
	perc_genes = 100.*perc_genes/float(len(gene2go))
	#print(perc_genes, float(len(gene2go)))    
	return [BH_enrichment, MNE, perc_genes, perc_enr_cluster, enriched_genes]

def GO_meaning(clusters_GOs):
    GOterms_meaning = []
    for clust in clusters_GOs:
        GOterms_meaning_clst = []
        for GO_term in clust:
            try:
                GOterm_meaning = godag[GO_term].name
            except KeyError:
                GOterm_meaning = GO_term
            GOterms_meaning_clst.append(GOterm_meaning)
        GOterms_meaning.append(GOterms_meaning_clst)    
    return GOterms_meaning

def parseKEGGnames(KEGGnames_file):
    KEGG2names = {}
    for line in KEGGnames_file.readlines():
        lspt=line.strip().split('\t')
        KEGG_name = lspt[0]
        name = lspt[1]
        KEGG2names[KEGG_name] = name
    return(KEGG2names)    

def ReactomeNames(ReactomeNamesfile):
    Reactome2names = {}
    for line in ReactomeNamesfile.readlines():
        lspt=line.strip().split('\t')
        Reactome_name = lspt[0]
        name = lspt[1]
        Reactome2names[Reactome_name] = name
    return(Reactome2names) 

def write_txt_enrRes(enrRes_dict, Clusters, sd, save_file, save_ext):
    with open(f'{sd}/{save_file}_{save_ext}.txt', 'w') as f: 
        for i in enrRes_dict:
            f.write(f'{i}\t({len(Clusters[i])} genes)\n')
            for pair in enrRes_dict[i]:                        
                f.write(f'{pair[0]}\t{pair[1]}\n')
            f.write('\n')
 
def EnrRes2Namesdict(EnrRes, term2name):
    EnrClsts_dict = {i:[] for i in range(len(EnrRes))}
    for i in range(len(EnrRes)):
        clst = EnrRes[i]
        for pair in clst:
            name = term2name[pair[0]]
            EnrClsts_dict[i].append([name,pair[1]])       
    return EnrClsts_dict
    
            
def Sort(sub_li):
 
    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of
    # sublist lambda has been used
    sub_li.sort(key = lambda x: x[1])
    return sub_li

def plot_nterms_across_clusters(terms_enrichment_clsts, n_clusters_range, title, save_file):
    terms_enrichment_clsts_l = {k:[] for k in n_clusters_range}   
    for k1 in terms_enrichment_clsts:
        for k2 in terms_enrichment_clsts[k1]:
            terms_enrichment_clsts_l[k1].append(len(terms_enrichment_clsts[k1][k2]))

    fig, ax = plt.subplots(figsize=(8, 5))
    for k1 in terms_enrichment_clsts_l:
        for l in terms_enrichment_clsts_l[k1]:        
            plt.scatter(k1,l)
    ax.set_title(title)
    ax.set_ylabel('Number of Enriched terms')
    ax.set_xticks(range(1,len(n_clusters_range)+1), n_clusters_range)
    plt.savefig(save_file,  dpi = 350, format='jpg', bbox_inches='tight')	
    plt.show()
    plt.close()
    
"""
	Main code starts here
"""


KEGGnames_file = open('input/HSA_KEGG_Pathways_meaning.lst', 'r')
KEGG2names = parseKEGGnames(KEGGnames_file)
Reactomefilenames = open('input/HSA_Reactome_Pathways_meaning.lst', 'r')
Reactome2names = ReactomeNames(Reactomefilenames)
ReactomeGroups = pd.read_csv('input/Reactome_Group.csv')
ReactomeRoot_dict = dict(zip(ReactomeGroups['Pathway'], ReactomeGroups['Group']))
ReactomeSubgroups_dict = dict(zip(ReactomeGroups['Pathway'], ReactomeGroups['Subgroup']))
ReactomeSubgroups2Root_dict = dict(zip(ReactomeGroups['Subgroup'], ReactomeGroups['Group']))


All_genes_inters = []
for file in os.listdir('input/Genelists'):
    Geneslist = pd.read_csv(f'input/Genelists/{file}')
    Geneslist = Geneslist['genes'].tolist()
    All_genes_inters.append(Geneslist)
        
All_genes_inters = list(set.intersection(*map(set,All_genes_inters)))
fgoname = 'input/PDmap_BasePathways_noSBUK.lst'
PDmap_counts, gene2PDmap = load_PD(fgoname, All_genes_inters)
print(f'{len(gene2PDmap)} PDmap_BasePathways_noSBUK genes')

fgoname_kp = 'input/HSA_Kegg_Pathways.lst'
go_counts_kp,gene2go_kp = load_GO(fgoname_kp, All_genes_inters)
print("%i KP annotated genes"%(len(gene2go_kp)))
   
fgoname_rp = 'input/HSA_Reactome_Pathways.lst'
go_counts_rp,gene2go_rp = load_GO(fgoname_rp, All_genes_inters)
print("%i RP annotated genes"%(len(gene2go_rp)))
   
fgoname_bp = 'input/HSA_GO-BP.lst'
go_counts_bp,gene2go_bp = load_GO(fgoname_bp, All_genes_inters)
print("%i BP annotated genes"%(len(gene2go_bp)))    
    
Go2name = {}
for goterm in  go_counts_bp:
    Go2name[goterm] = godag[goterm].name
    
n_clusters_range = list(range(1,21))

   
with open('output/PINK1_1neighCorePDPreds_PPI.txt') as f:    
    PPI_firstNeigh = f.read().splitlines()[:-1] 
Clusters = [PPI_firstNeigh]
sd = 'output'
save_file = 'PINK1_1neighCorePDPreds_PPI'

bh_PDmap, mne_PDmap, eg_PDmap, perc_cluster_PDmap, enriched_genes_PDmap = go_enrichment(Clusters, gene2PDmap, PDmap_counts)
# print(bh_PDmap)
bh_kp, mne_kp, eg_kp, perc_cluster_kp, enriched_genes_kp = go_enrichment(Clusters, gene2go_kp, go_counts_kp)            
bh_rp, mne_rp, eg_rp, perc_cluster_rp, enriched_genes_rp = go_enrichment(Clusters, gene2go_rp, go_counts_rp)
bh_bp, mne_bp, eg_bp, perc_cluster_bp, enriched_genes_bp = go_enrichment(Clusters, gene2go_bp, go_counts_bp)        


# write PDmap
EnrClsts_PDmap = {v: Sort(k) for v, k in enumerate(bh_PDmap)}
write_txt_enrRes(EnrClsts_PDmap, Clusters, sd, save_file, save_ext='PDmap')
 

EnrClsts_Kegg = {v: Sort(k) for v, k in enumerate(bh_kp)}
write_txt_enrRes(EnrClsts_Kegg, Clusters, sd, save_file, save_ext='KPterms')
EnrClsts_KeggNames = EnrRes2Namesdict(bh_kp, KEGG2names)        
write_txt_enrRes(EnrClsts_KeggNames, Clusters, sd, save_file, save_ext='KPmeaning')
# KP_enrichment_clsts[n_clusters] = EnrClsts_KeggNames
      

#write Reactome terms
EnrClsts_RP= {v: Sort(k) for v, k in enumerate(bh_rp)}
write_txt_enrRes(EnrClsts_RP, Clusters, sd, save_file, save_ext='RPterms')
EnrClsts_RPNames = EnrRes2Namesdict(bh_rp, Reactome2names)        
write_txt_enrRes(EnrClsts_RPNames, Clusters, sd, save_file, save_ext='RPmeaning')
# RP_enrichment_clsts[n_clusters] = EnrClsts_RPNames


#write GO terms
EnrClsts_BP= {v: Sort(k) for v, k in enumerate(bh_bp)}
write_txt_enrRes(EnrClsts_BP, Clusters, sd, save_file, save_ext='BPterms')
EnrClsts_BPNames = EnrRes2Namesdict(bh_bp, Go2name)        
write_txt_enrRes(EnrClsts_BPNames, Clusters, sd, save_file, save_ext='BPmeaning')
# BP_enrichment_clsts[n_clusters] = EnrClsts_BPNames


   




