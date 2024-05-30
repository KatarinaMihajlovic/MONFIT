# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 14:01:07 2023

@author: kmihajlo
"""
import math, json
from scipy.stats import hypergeom
import numpy as np
import matplotlib.pyplot as plt
import os.path
import pandas as pd
from matplotlib_venn import venn2, venn2_circles
from matplotlib.cm import ScalarMappable
import pickle

from goatools.obo_parser import GODag
godag = GODag(obo_file="input/go-basic.obo") 

with open('input/Entrez2Sym.pkl', 'rb') as handle:
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


def load_PDGO(fname, genes):
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
            try:
                name = term2name[pair[0]]
            except Exception:
                name = pair[0]
            EnrClsts_dict[i].append([name,pair[1]])       
    return EnrClsts_dict
    
            
def Sort(sub_li):
 
    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of
    # sublist lambda has been used
    sub_li.sort(key = lambda x: x[1])
    return sub_li


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

perc = '1.094'

GM_MO_df = pd.read_csv('input/MO_GM_D8_D18_D25_D32_D37_1.533_PDpreds.csv', header = 0, index_col=0)
GM_NB_df = pd.read_csv(f'input/NB_GM_D8_D18_D37_{perc}_PDpreds.csv', header = 0, index_col=0)

MO_preds = list(GM_MO_df.index)
NB_preds = list(GM_NB_df.index)
        
inters = list(set(NB_preds) & set(MO_preds))
Preds_NB_only = list(set(NB_preds) - set(MO_preds)) 
Preds_MO_only = list(set(MO_preds) - set(NB_preds)) 

print(len(inters))
print(len(Preds_MO_only))
print(len(Preds_NB_only))

TGM_MO_df = pd.read_csv('input/MO_GM_Length_D8_D18_D25_D32_D37_1.533_PDpreds.csv', header = 0, index_col=0)
TGM_MO_df['Rank'] =  np.arange(1, len(TGM_MO_df) + 1)
TGM_MO_only = TGM_MO_df.loc[Preds_MO_only].sort_values(by='Rank')


days = [8, 18, 37]
All_genes_inters = []
for file in os.listdir('input/Genelists'):
    day = file.split('_')[2][:-4]
    if int(day) in days:
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
    
# fgoname = 'input/PD_GO_UCL.lst'
# PDgo_counts_bp,PDgene2go_bp = load_PDGO(fgoname, All_genes_inters)
# PD_GOs = set()
# fRead = open(fgoname, "r")
# for line in fRead.readlines():
#     lspt = line.strip().split()
#     go = lspt[-1]
#     PD_GOs.add(go)


Go2name = {}
for goterm in  go_counts_bp:
    Go2name[goterm] = godag[goterm].name


genes_lists = {'Preds_NB':NB_preds, 'Preds_MO':MO_preds, 'inters':inters, 'Preds_NB_only':Preds_NB_only, 'Preds_MO_only':Preds_MO_only}
# genes_lists = {'Preds_MO':Preds_MO}
for save_file in genes_lists:
    print(save_file)
    Clusters = [genes_lists[save_file]]
    sd = f'output/{save_file}'    
    if not os.path.exists(sd):
        os.makedirs(sd)
          
    
    bh_PDmap, mne_PDmap, eg_PDmap, perc_cluster_PDmap, enriched_genes_PDmap = go_enrichment(Clusters, gene2PDmap, PDmap_counts)
    # bh_PDGO, mne_PDGO, eg_PDGO, perc_cluster_PDGO, enriched_genes_PDGO = go_enrichment(Clusters, PDgene2go_bp, PDgo_counts_bp)

    # print(bh_PDmap)
    bh_kp, mne_kp, eg_kp, perc_cluster_kp, enriched_genes_kp = go_enrichment(Clusters, gene2go_kp, go_counts_kp)            
    bh_rp, mne_rp, eg_rp, perc_cluster_rp, enriched_genes_rp = go_enrichment(Clusters, gene2go_rp, go_counts_rp)
    bh_bp, mne_bp, eg_bp, perc_cluster_bp, enriched_genes_bp = go_enrichment(Clusters, gene2go_bp, go_counts_bp)        
    
    
    # write PDmap
    EnrClsts_PDmap = {v: Sort(k) for v, k in enumerate(bh_PDmap)}
    write_txt_enrRes(EnrClsts_PDmap, Clusters, sd, save_file, save_ext='PDmap')
 
    # # write PD GO
    # EnrClsts_PDGO= {v: Sort(k) for v, k in enumerate(bh_PDGO)}
    # write_txt_enrRes(EnrClsts_PDGO, Clusters, sd, save_file, save_ext='PDGO')
    # EnrClsts_PDGONames = EnrRes2Namesdict(bh_PDGO, Go2name)        
    # write_txt_enrRes(EnrClsts_PDGONames, Clusters, sd, save_file, save_ext='PDGOmeaning') 
    
    # write KEGG
    EnrClsts_Kegg = {v: Sort(k) for v, k in enumerate(bh_kp)}
    write_txt_enrRes(EnrClsts_Kegg, Clusters, sd, save_file, save_ext='KPterms')
    EnrClsts_KeggNames = EnrRes2Namesdict(bh_kp, KEGG2names)        
    write_txt_enrRes(EnrClsts_KeggNames, Clusters, sd, save_file, save_ext='KPmeaning')
           
    
    #write Reactome terms
    EnrClsts_RP= {v: Sort(k) for v, k in enumerate(bh_rp)}
    write_txt_enrRes(EnrClsts_RP, Clusters, sd, save_file, save_ext='RPterms')
    EnrClsts_RPNames = EnrRes2Namesdict(bh_rp, Reactome2names)        
    write_txt_enrRes(EnrClsts_RPNames, Clusters, sd, save_file, save_ext='RPmeaning')
    
    
    #write GO terms
    EnrClsts_BP= {v: Sort(k) for v, k in enumerate(bh_bp)}
    write_txt_enrRes(EnrClsts_BP, Clusters, sd, save_file, save_ext='BPterms')
    EnrClsts_BPNames = EnrRes2Namesdict(bh_bp, Go2name)        
    write_txt_enrRes(EnrClsts_BPNames, Clusters, sd, save_file, save_ext='BPmeaning')
    #PD GO terms
    enrGOs = bh_bp[0]
    enrGOs = [x[0] for x in enrGOs]
    # enrPDgos = list(set(enrGOs) & set(PD_GOs))
    # enrPDgos_names = [Go2name[x] for x in enrPDgos]
    # print(len(enrPDgos))
    # print(enrPDgos_names)

    # save PD preds that are annotated with enriched Kegg pathways; save Pathways and Core PD preds that annotate them 
    Clusters_flatten = [item for sublist in Clusters for item in sublist]
    KP_terms = []
    KP_terms_meaning_p = []
    for clust in bh_kp:
        for pair in clust:
            try:
                term = pair[0]
                KP_terms.append(term)
                KP_terms_meaning_p.append([term, pair[1]])
            except IndexError:
                pass
                
    PD_preds_KEGGEnr = []
    Unannotated_genes = []
    for gene in Clusters_flatten:
        try:
            terms= gene2go_kp[gene]
            for term in terms:
                if term in KP_terms:
                    PD_preds_KEGGEnr.append(gene)
        except KeyError:
            Unannotated_genes.append(gene)

    PD_preds_KEGGEnr = list(set(PD_preds_KEGGEnr))
    with open(f'{sd}/{save_file}_KPEnirchedGenes.txt', 'w') as f: 
        for gene in PD_preds_KEGGEnr:
            f.write(f'{gene}\n')   
            

    with open(f'{sd}/{save_file}_KPEnirchedTerms_Genes.txt', 'w') as f: 
        terms_genes = []
        terms_genes_count_p = []
        for term in KP_terms:
            term_genes = []
            # print(term, KEGG2names[term])
            # f.write(f'{KEGG2names[term]}\n')   
            for gene, KP_terms_d in gene2go_kp.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
                if term in KP_terms_d and gene in Clusters_flatten:
                    term_genes.append(gene)                                # f.write(f'{gene}\n')   
            terms_genes.append([term, term_genes])
            
            for term_p in KP_terms_meaning_p:
                # print(term_p, term)
                if KEGG2names[term] == term_p[0]:
                    terms_genes_count_p.append([KEGG2names[term][:-3], term_p[1], len(term_genes)])
        
        terms_genes = Sort(terms_genes)
        for i in range(len(terms_genes)):
            f.write(f'{KEGG2names[terms_genes[i][0]]}\t{terms_genes[i][0]}\t{len(terms_genes[i][1])}\n')
            for j in range(len(terms_genes[i][1])):
                f.write(f'{terms_genes[i][1][j]}\n')
            f.write('\n')
            
    #reactome Subgroups  
    if len(Clusters) == 1:
        TopPreds_RP_Root = {key:[] for key in set(ReactomeGroups['Group'])}
        TopPreds_RP_Subgroup = {key:[] for key in set(ReactomeGroups['Subgroup'])}
        reactomes_enr = bh_rp[0]
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
        
        with open(f'{sd}/{save_file}_ReactomeRoots.txt', 'w') as f:
            for root_path in TopPreds_RP_Root:
                f.write(f'{root_path}\t{len(TopPreds_RP_Root[root_path])}\n')
                for path in TopPreds_RP_Root[root_path]:
                    f.write(f'{path}\n')
                f.write('\n')
        
        Root2SG_Preds_size = [[root,len(TopPreds_RP_Root[root])] for root in TopPreds_RP_Root]
        Sort(Root2SG_Preds_size)
        order  = [x[0] for x in Root2SG_Preds_size]
        
        with open(f'{sd}/{save_file}_ReactomeRootsSubGroups.txt', 'w') as f:
            for root_path in order:
                f.write(f'{root_path}\t{len(Root2SG_Preds[root_path])}\t({len(TopPreds_RP_Root[root_path])})\n')
                for subgroup_d in Root2SG_Preds[root_path]:
                    # print(root_path, subgroup_d)
                    for subgroup in subgroup_d:
                        f.write(f'\t{subgroup}\t{len(subgroup_d[subgroup])}\n')    
                        for path in subgroup_d[subgroup]:
                            f.write(f'\t\t{path}\n')
                f.write('\n')      
    
  
    
import venn
# compre DEGS, DAPS, MONFIT enriched pathways
EPs = ['RP', 'BP'] 
EP_D = {'RP':'RP', 'BP':'GO-BP'}

in_dir = 'output'

for EP in EPs:
    print(EP)
    with open(f'{in_dir}/Preds_MO_only/Preds_MO_only_{EP}terms.txt') as f:
        f.readline()
        EP_MONFIT = f.read().splitlines() 
    EP_MONFIT = [x.split('\t')[0] for x in EP_MONFIT][:-1]    
 
    # with open(f'{in_dir}/Preds_MO_only/Pred_MO_only_{EP}meaning.txt') as f:
    #     f.readline()
    #     EP_MONFIT_meaning = f.read().splitlines() 
    # EP_MONFIT_meaning = [x.split('\t')[0] for x in EP_MONFIT_meaning][:-1]    
    

        
    with open(f'{in_dir}/Preds_NB_only/Preds_NB_only_{EP}terms.txt') as f:
        f.readline()
        EP_NBs = f.read().splitlines() 
    EP_NBs = [x.split('\t')[0] for x in EP_NBs][:-1]    
        
    print(list(set(EP_NBs)&set(EP_MONFIT)))    

    # Venn
    genes = [set(EP_MONFIT), set(EP_NBs)]
    gene_sets = ['MO', 'NB']
    labels = venn.generate_petal_labels(genes)#, fill=['number', 'logic'])
    fig, ax = venn.venn2(labels, names=gene_sets, fontsize =24, figsize= (9, 12))
    ax.set_title(EP_D[EP], fontsize = 28, y=0.9)  
    if EP == 'BP':
        leg = ax.legend(gene_sets, bbox_to_anchor=(0.9, 1), fancybox=True, fontsize = 22)
    else:
        ax.get_legend().remove()

    plt.savefig(f'output/{EP}_Overlap.jpg',  dpi = 350, format='jpg', bbox_inches='tight')	
    plt.show()    
  

  

  

  

  
    
  
    
  


