# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:15:35 2023

@author: kimhajlo
"""

#import networkx as nx
import math, sys
from scipy.stats import hypergeom
import numpy as np
import os.path
import pandas as pd
import pickle
from goatools.obo_parser import GODag

wd = str(sys.argv[1])
godag = GODag(obo_file=f'{wd}/input/go-basic.obo') 

with open(f'{wd}/input/Entrez2Sym.pkl', 'rb') as handle:
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
  

def ReactomeNames(ReactomeNamesfile):
    Reactome2names = {}
    for line in ReactomeNamesfile.readlines():
        lspt=line.strip().split('\t')
        Reactome_name = lspt[0]
        name = lspt[1]
        Reactome2names[Reactome_name] = name
    ReactomeNamesfile.close()
    return(Reactome2names) 

def ReactomeNames_inv(ReactomeNamesfile):
    Reactome2names = {}
    for line in ReactomeNamesfile.readlines():
        lspt=line.strip().split('\t')
        Reactome_name = lspt[0]
        name = lspt[1]
        Reactome2names[name] = Reactome_name
    ReactomeNamesfile.close()
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

 
def GenesEnrichedTerms(Clusters, enriched_terms, annot_name, annot, annot2names, enriched_genes):     #   annot name gene2go_kp                 
    Clusters_flatten = [item for sublist in Clusters for item in sublist]
    KP_terms = []
    KP_terms_meaning_p = []
    for clust in enriched_terms:
        for pair in clust:
            try:
                term = pair[0]
                KP_terms.append(term)
                KP_terms_meaning_p.append([term, pair[1]])
            except IndexError:
                pass
                

    try:
        enriched_genes = [item for sublist in enriched_genes for item in sublist]
    except Exception:
        pass

    PD_preds_KEGGEnr = list(set(enriched_genes))
    with open(f'{sd}/{save_file}_{annot}EnirchedGenes.txt', 'w') as f: 
        for gene in PD_preds_KEGGEnr:
            f.write(f'{gene}\n')   
            

    with open(f'{sd}/{save_file}_{annot}EnirchedTerms_Genes.txt', 'w') as f: 
        terms_genes = []
        terms_genes_count_p = []
        for term in KP_terms:
            term_genes = []
            # print(term, KEGG2names[term])
            # f.write(f'{KEGG2names[term]}\n')   
            for gene, KP_terms_d in annot_name.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
                if term in KP_terms_d and gene in Clusters_flatten:
                    term_genes.append(gene)                                # f.write(f'{gene}\n')   
            terms_genes.append([term, term_genes])
            
            for term_p in KP_terms_meaning_p:
                # print(term_p, term)
                if annot2names[term] == term_p[0]:
                    terms_genes_count_p.append([annot2names[term][:-3], term_p[1], len(term_genes)])
        
        terms_genes = Sort(terms_genes)
        for i in range(len(terms_genes)):
            f.write(f'{annot2names[terms_genes[i][0]]}\t{terms_genes[i][0]}\t{len(terms_genes[i][1])}\n')
            for j in range(len(terms_genes[i][1])):
                f.write(f'{terms_genes[i][1][j]}\n')
            f.write('\n')   

"""
	Main code starts here
"""


Reactomefilenames = open(f'{wd}/input/HSA_Reactome_Pathways_meaning.lst', 'r')
Reactome2names = ReactomeNames(Reactomefilenames)
Reactomefilenames = open(f'{wd}/input/HSA_Reactome_Pathways_meaning.lst', 'r')
Names2RP =  ReactomeNames_inv(Reactomefilenames)
ReactomeGroups = pd.read_csv(f'{wd}/input/Reactome_Group.csv')
ReactomeRoot_dict = dict(zip(ReactomeGroups['Pathway'], ReactomeGroups['Group']))
ReactomeSubgroups_dict = dict(zip(ReactomeGroups['Pathway'], ReactomeGroups['Subgroup']))
ReactomeSubgroups2Root_dict = dict(zip(ReactomeGroups['Subgroup'], ReactomeGroups['Group']))

unique_genes = []

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
       
    fgoname_rp = f'{wd}/input/HSA_Reactome_Pathways.lst'
    go_counts_rp,gene2go_rp = load_GO(fgoname_rp, All_genes_inters)
    print("%i RP annotated genes"%(len(gene2go_rp)))
       
    fgoname_bp = f'{wd}/input/HSA_GO-BP.lst'
    go_counts_bp,gene2go_bp = load_GO(fgoname_bp, All_genes_inters)
    print("%i BP annotated genes"%(len(gene2go_bp)))    
        
    Go2name = {}
    for goterm in  go_counts_bp:
        Go2name[goterm] = godag[goterm].name
        

    
    for perc_case in os.listdir(f'{in_dir}/{days_set}'):     
        for file in os.listdir(f'{in_dir}/{days_set}/{perc_case}'):
                            
            save_file = file[:-4]
            print(file)
            sd = f'{wd}/output/{days_set}/{perc_case}'
            if not os.path.exists(sd):
                os.makedirs(sd)            
         
            if file.endswith('pkl'):
                with open(f'{in_dir}/{days_set}/{perc_case}/{file}', 'rb') as handle:
                    Clusters = pickle.load(handle)
            elif file.endswith('txt'):
                with open(f'{in_dir}/{days_set}/{perc_case}/{file}') as f:
                    lines = f.readlines()
                Clusters = [[line[:-1] for line in lines]]
            elif file.endswith('csv'):
                Clusters = pd.read_csv(f'{in_dir}/{days_set}/{perc_case}/{file}', index_col=0)
                Clusters = [list(Clusters.index)]
                
            bh_rp, mne_rp, eg_rp, perc_cluster_rp, enriched_genes_rp = go_enrichment(Clusters, gene2go_rp, go_counts_rp)
            bh_bp, mne_bp, eg_bp, perc_cluster_bp, enriched_genes_bp = go_enrichment(Clusters, gene2go_bp, go_counts_bp)        
            
                               
            
            #write Reactome terms
            EnrClsts_RP= {v: Sort(k) for v, k in enumerate(bh_rp)}
            write_txt_enrRes(EnrClsts_RP, Clusters, sd, save_file, save_ext='RPterms')
            EnrClsts_RPNames = EnrRes2Namesdict(bh_rp, Reactome2names)        
            write_txt_enrRes(EnrClsts_RPNames, Clusters, sd, save_file, save_ext='RPmeaning')
            # EnrClsts_RP_REP  = EnrClsts_RP[0]
            
                  
            
            #write GO terms
            EnrClsts_BP= {v: Sort(k) for v, k in enumerate(bh_bp)}
            write_txt_enrRes(EnrClsts_BP, Clusters, sd, save_file, save_ext='BPterms')
            EnrClsts_BPNames = EnrRes2Namesdict(bh_bp, Go2name)        
            write_txt_enrRes(EnrClsts_BPNames, Clusters, sd, save_file, save_ext='BPmeaning')
            
        
            # # save Reactome and GO enriched terms per cluster in individual files
            sdRP = f'{sd}/ClustsEnrichRP'
            if not os.path.exists(sdRP):
                os.makedirs(sdRP)      
            sdBP = f'{sd}/ClustsEnrichBP'
            if not os.path.exists(sdBP):
                os.makedirs(sdBP)     
            
            for k in EnrClsts_RP:
                with open(f'{sdRP}/RPterms_{file}.txt','w') as f:
                    for term in EnrClsts_RP[k]:
                        f.write(f'{term[0]}\n')
            for k in EnrClsts_BP:
                with open(f'{sdBP}/BPterms_{file}.txt','w') as f:
                    for term in EnrClsts_BP[k]:
                        f.write(f'{term[0]}\n')                       
       
        
            GenesEnrichedTerms(Clusters, bh_rp, gene2go_rp, 'RP', Reactome2names, enriched_genes_rp)
            GenesEnrichedTerms(Clusters, bh_bp, gene2go_bp, 'BP', Go2name, enriched_genes_bp)

            
            Clusters_flatten = [item for sublist in Clusters for item in sublist]
            
            #reactome Subgroups  
            # if len(Clusters) == 1:
            sdRPSub = f'{sd}/ReactomeSubgroups'
            if not os.path.exists(sdRPSub):
                os.makedirs(sdRPSub)                           

            for i in range(len(bh_rp)):
                TopPreds_RP_Root = {key:[] for key in set(ReactomeGroups['Group'])}
                TopPreds_RP_Subgroup = {key:[] for key in set(ReactomeGroups['Subgroup'])}
                reactomes_enr = bh_rp[i]
                reactomes_enr = [x[0] for x in reactomes_enr]
                for reactome in reactomes_enr:
                    TopPreds_RP_Root[ReactomeRoot_dict[reactome]].append(Reactome2names[reactome])
                    TopPreds_RP_Subgroup[ReactomeSubgroups_dict[reactome]].append(Reactome2names[reactome])
                
                TopPreds_RP_Root = {k:v for k,v in TopPreds_RP_Root.items() if v}
                TopPreds_RP_Subgroup = {k:v for k,v in TopPreds_RP_Subgroup.items() if v}
                
                TopPreds_RP_Root_summary = {k:v for k,v in TopPreds_RP_Root.items() if len(v)>1}
                TopPreds_RP_Subgroup_summary = {k:v for k,v in TopPreds_RP_Subgroup.items() if len(v)>1}
                
                Root2SG_Preds = {key:[] for key in set(TopPreds_RP_Root.keys())}
                for SG in TopPreds_RP_Subgroup:
                    d = {}
                    d[SG] = TopPreds_RP_Subgroup[SG]
                    Root2SG_Preds[ReactomeSubgroups2Root_dict[SG]].append(d)
                
                with open(f'{sdRPSub}/GM_Length_{days_set}_{perc_case}_PDpreds_{i}_ReactomeRoots.txt', 'w') as f:
                    for root_path in TopPreds_RP_Root:
                        f.write(f'{root_path}\t{len(TopPreds_RP_Root[root_path])}\n')
                        for path in TopPreds_RP_Root[root_path]:
                            f.write(f'{path}\n')
                        f.write('\n')
                
                Root2SG_Preds_size = [[root,len(TopPreds_RP_Root[root])] for root in TopPreds_RP_Root]
                Sort(Root2SG_Preds_size)
                order  = [x[0] for x in Root2SG_Preds_size]
                
                with open(f'{sdRPSub}/GM_Length_{days_set}_{perc_case}_PDpreds_{i}_ReactomeRootsSubGroups.txt', 'w') as f:
                    for root_path in order:
                        f.write(f'{root_path}\t{len(Root2SG_Preds[root_path])}\t({len(TopPreds_RP_Root[root_path])})\n')
                        for subgroup_d in Root2SG_Preds[root_path]:
                            # print(root_path, subgroup_d)
                            for subgroup in subgroup_d:
                                f.write(f'\t{subgroup}\t{len(subgroup_d[subgroup])}\n')    
                                for path in subgroup_d[subgroup]:
                                    f.write(f'\t\t{path}\n')
                        f.write('\n')       
                        
    
   