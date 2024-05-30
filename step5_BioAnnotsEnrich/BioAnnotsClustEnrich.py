import math, sys, pickle
from scipy.stats import hypergeom
import numpy as np
import os.path
import pandas as pd
# import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag



E_w = float(sys.argv[1])
PPI_w = float(sys.argv[2])
COEX_w = float(sys.argv[3])
GI_w = float(sys.argv[4])
MI_w = float(sys.argv[5])
wd = str(sys.argv[6])


# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)


with open(f'{wd}/input/Entrez2Sym.pkl', 'rb') as handle:
    Entrez2Sym = pickle.load(handle)
    
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p


# Read GO terms annotations
def load_GO(fname, genes):
    geneset=set(genes)
    GO_counter = {}
    gene2GO = {}
    #loading raw gene & go relationships
    fRead = open(fname, "r")
    for line in fRead.readlines():
        lspt = line.strip().split()
        # print(lspt)
        if len(lspt)>1:
            try:
                geneid = Entrez2Sym[lspt[0]]
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



def avg_percentile_enrich(ER_ccs):
    for i in range(len(ER_ccs)):
        mean_val = np.percentile(ER_ccs[i], 50)
        ER_ccs[i].append(mean_val)
        
        q3, q1 = np.percentile(ER_ccs[i], [84.13, 15.87])
        min_error = mean_val - q1
        max_error = q3 - mean_val 
        ER_ccs[i].append(min_error)
        ER_ccs[i].append(max_error)

    return ER_ccs


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return zip(*[iter(iterable)]*n)

def plotGeneClusters(tick_labels, ER_list_BP_ccs, ER_list_KP_ccs, ER_list_RP_ccs, clst_type, comb, save_dir):
    capsize = 5 
    alpha = 0.7
    N = len(tick_labels)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2       # the width of the bars
   
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    for i in range(N):
        # print(all_nets[i-1], ae1, ER_list_KP_ccs[i-1])
        rects1 = ax.bar(ind[i-1], ER_list_KP_ccs[i-1], width=width, color='r', alpha = alpha, ecolor='black', capsize=capsize)
        rects2 = ax.bar(ind[i-1]+width, ER_list_RP_ccs[i-1], width=width, color='g', alpha = alpha, ecolor='black', capsize=capsize)
        rects3 = ax.bar(ind[i-1]+width*2, ER_list_BP_ccs[i-1], width=width, color='b', alpha = alpha, ecolor='black', capsize=capsize)

        # print(all_nets[i-1])
        # print(ER_clusters_KP[i-1])
        # print(ER_clusters_RP[i-1])
        # print(ER_clusters_BP[i-1])
        # print('\n')    
        
    ax.set_ylabel(f'{clst_type} with enriched annotations (%)', fontsize = 24, fontweight = 'bold')
    ax.set_xticks(ind+width)
    # if labels_key[cell_cond] == 'C21' or labels_key[cell_cond] == 'PD21':
    tick_labels_sup = []
    for tick in tick_labels:
        splt = tick.split('_')
        tick_sup = f'{splt[0]}$_{{{splt[1]}}}$'
        tick_labels_sup.append(tick_sup)
        
    ax.set_xticklabels(tick_labels_sup,  fontsize = 26) #nets_everythinga should be all_nets
    # else:
    #     ax.set_xticklabels(len(all_nets)*[''])
    title_sup = ''
    comb_splt = comb.split('_')
    for x, y in grouped(comb_splt, 2):
        pair =  f'{x}$_{{{y}}}$'
        title_sup+=pair
    
    plt.title(title_sup,  fontsize = 34, pad = 24)


    ax.tick_params(axis='y', which='major', labelsize=24)
    ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO-BP'),  fontsize = 24)

    if clst_type == 'Clusters':       
        plt.ylim(0,130)
        plt.tight_layout()    
        plt.savefig(f'{save_dir}/{comb}_ECs.jpg', dpi = 350, format='jpg')	
    elif clst_type == 'Genes':   
        plt.ylim(0,65)
        plt.tight_layout()    
        plt.savefig(f'{save_dir}/{comb}_EGs.jpg', dpi = 350, format='jpg')	
    plt.show()	
    plt.close()
    
def save_pkl(variable, sf, sd):
    if not os.path.exists(sd):
        os.makedirs(sd)  
    with open(f'{sd}/{sf}', 'wb') as handle:
        pickle.dump(variable, handle)   
 



       	        
"""
	Main code starts here
"""


godag = GODag(obo_file=f'{wd}/input/go-basic.obo') 
in_dir = 'step3_ClusterG1NMTF_WeightedNets/output'
comb = f'P_{PPI_w}_C_{COEX_w}_G_{GI_w}_M_{MI_w}_E_{E_w}'



Flag = True
if Flag == True:
    for cc in os.listdir(in_dir):
        print(cc)
        sd = f'{wd}/output/{cc}/{comb}'
        try:
            os.makedirs(sd)
        except FileExistsError:
            # directory already exists
            pass
        
        genes = pd.read_csv(f'{wd}/input/Genelists/Geneslist_{cc}.csv', header = 0)
        PPIlist = genes['genes'].to_list()
    
        fgoname_kp = f'{wd}/input/HSA_Kegg_Pathways.lst'
        go_counts_kp,gene2go_kp = load_GO(fgoname_kp, PPIlist)
        print(f'{len(gene2go_kp)} KP annotated genes')
        
        fgoname_rp = f'{wd}/input/HSA_Reactome_Pathways.lst'
        go_counts_rp,gene2go_rp = load_GO(fgoname_rp, PPIlist)
        print(f'{len(gene2go_rp)} RP annotated genes')
        
        fgoname_bp = f'{wd}/input/HSA_GO-BP.lst'
        go_counts_bp,gene2go_bp = load_GO(fgoname_bp, PPIlist)
        print(f'{len(gene2go_bp)} BP annotated genes')
        
        ER_clusters_KP = []
        ER_clusters_RP = []
        ER_clusters_BP = []
        
        ER_genes_KP = []
        ER_genes_RP = []
        ER_genes_BP = []   
    
        
        runs = []
        for file in os.listdir(f'{in_dir}/{cc}/{comb}'):
            print(file)
            run = file.split('_')[0]
            runs.append(run)
    
            with open(f'{in_dir}/{cc}/{comb}/{file}', 'rb') as handle:
                G1clusters = pickle.load(handle)
            Cp = [[str(y) for y in x] for x in G1clusters]           
           
            bh_kp, mne_kp, eg_kp, perc_cluster_kp = go_enrichment(Cp, gene2go_kp, go_counts_kp)
            bh_rp, mne_rp, eg_rp, perc_cluster_rp = go_enrichment(Cp, gene2go_rp, go_counts_rp)
            bh_bp, mne_bp, eg_bp, perc_cluster_bp = go_enrichment(Cp, gene2go_bp, go_counts_bp)
           					#print(f'Function Enirchment : {len(bh_bp)}')
           
            KP_terms = [[i[0] for i in nested] for nested in bh_kp]
            RP_terms = [[i[0] for i in nested] for nested in bh_rp]
            BP_terms = [[i[0] for i in nested] for nested in bh_bp]
           
            
            
            print(perc_cluster_kp, eg_kp)
            print(perc_cluster_rp, eg_rp)
            print(perc_cluster_bp, eg_bp)
           
            ER_clusters_KP.append(perc_cluster_kp)
            ER_clusters_RP.append(perc_cluster_rp)
            ER_clusters_BP.append(perc_cluster_bp)
           
            ER_genes_KP.append(eg_kp)
            ER_genes_RP.append(eg_rp)
            ER_genes_BP.append(eg_bp)
            
    
                                
        clst_enrich_df = pd.DataFrame(list(zip(ER_clusters_BP, ER_clusters_KP, ER_clusters_RP)), index = runs, columns =['BP', 'KP', 'RP'])
        gene_enrich_df = pd.DataFrame(list(zip(ER_genes_BP, ER_genes_KP, ER_genes_RP)), index = runs, columns =['BP', 'KP', 'RP'])    
        
        clst_enrich_df.to_csv(f'{sd}/EC.csv')
        gene_enrich_df.to_csv(f'{sd}/EG.csv')


index_ord = ['WT_8','WT_18','WT_25','WT_32','WT_37','ND_8','ND_18','ND_25','ND_32','ND_37']

df_mean_EC = pd.DataFrame(columns=['BP', 'KP', 'RP'])
df_stdv_EC = pd.DataFrame(columns=['BP', 'KP', 'RP'])
df_mean_EG = pd.DataFrame(columns=['BP', 'KP', 'RP'])
df_stdv_EG = pd.DataFrame(columns=['BP', 'KP', 'RP'])


for cc in os.listdir(in_dir):
    print(cc)   
    sd = f'{wd}/output/{cc}'

    EC_df = pd.read_csv(f'{sd}/{comb}/EC.csv', index_col=0)
    mean = EC_df.mean()
    stdv = EC_df.std()
    df_mean_EC.loc[cc] = mean
    df_stdv_EC.loc[cc] = stdv
    
    
    EG_df = pd.read_csv(f'{sd}/{comb}/EG.csv', index_col=0)
    mean = EG_df.mean()
    stdv = EG_df.std()
    df_mean_EG.loc[cc] = mean
    df_stdv_EG.loc[cc] = stdv
 
df_mean_EC = df_mean_EC.reindex(index_ord)
df_mean_EG = df_mean_EG.reindex(index_ord)
   

# plot clusters
ax = df_mean_EC.plot(kind='bar', yerr=df_stdv_EC, capsize=5 ,figsize=(12, 8))
ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)

ax.set_xlabel('Cell Conditions', fontsize=20)
ax.set_ylabel('Clusters with enriched annotations (%)', fontsize=20)
plt.tick_params(axis='x', labelsize=16, rotation=0)
plt.tick_params(axis='y', labelsize=16)
plt.ylim(0,130)
ax.legend(fontsize = 16)

plt.savefig(f'{wd}/output/ECs.jpg', dpi = 350, format='jpg')	
plt.show()
plt.close()

# plot genes
ax = df_mean_EG.plot(kind='bar', yerr=df_stdv_EG, capsize=5 ,figsize=(12, 8))
ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)

ax.set_xlabel('Cell Conditions', fontsize=20)
ax.set_ylabel('Genes with enriched annotations (%)', fontsize=20)
plt.tick_params(axis='x', labelsize=16, rotation=0)
plt.tick_params(axis='y', labelsize=16)
plt.ylim(0,60)
ax.legend(fontsize = 16)

plt.savefig(f'{wd}/output/EGs.jpg', dpi = 350, format='jpg')	
plt.show()
plt.close()
    

print(df_mean_EG.mean())
    
    
    
    
    
    
    
