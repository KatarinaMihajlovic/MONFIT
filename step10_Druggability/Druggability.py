# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:50:58 2023

@author: kmihajlo
"""

import pandas as pd
import os, pickle
from collections import Counter
import math, sys
import numpy as np
from scipy.stats import hypergeom


	
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

def Sort(sub_li):
 
    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of
    # sublist lambda has been used
    sub_li.sort(key = lambda x: x[1])
    return sub_li

'''MAIN'''
wd = str(sys.argv[1])
 
All_genes_inters = []
for file in os.listdir(f'{wd}/input/Genelists'):
    Geneslist = pd.read_csv(f'{wd}/input/Genelists/{file}')
    Geneslist = Geneslist['genes'].tolist()
    All_genes_inters.append(Geneslist)
All_genes_inters = list(set.intersection(*map(set,All_genes_inters)))
           
    
DTI = pd.read_csv(f'{wd}/input/DTI_invest_approv.csv')
DTI = DTI.drop_duplicates(keep='first')
genes = list(DTI['Gene Target'].values)
drugs = list(DTI['DrugBank ID'].values)
drug_names = list(DTI['DrugBank Name'].values)
drug2name = dict(zip(drugs,drug_names))

print(len(set(drugs)))


dict2 = dict(zip(drugs,genes))
gene2drug = {}
drug_count = {}

for gene,drug in zip(genes,drugs):
    if gene not in gene2drug:
        gene2drug[gene] = [drug]
    else:
        gene2drug[gene].append(drug)
    if drug not in drug_count:
        drug_count[drug] = 1
    else:
        drug_count[drug] += 1


# Drugs_dict = {}

in_dir = f'{wd}/input/Predictions'
for days_set in os.listdir(in_dir): 
    # Drugs_dict[days_set] = {}  
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
    
    gene2drug_ds = gene2drug.copy()
    for gene in All_genes_inters:
        if gene not in gene2drug_ds:
            gene2drug_ds[gene] = []
        
      
    for perc_case in os.listdir(f'{in_dir}/{days_set}'):
        files = os.listdir(f'{in_dir}/{days_set}/{perc_case}')
        # files.append('GM_D18_D25_D32_D37.csv')
        for file in files:       
            # print(file)
        # perc_case    = '20perc'
        # Drugs_dict[perc_case] = {}      
            print(perc_case, file)
            save_file = file[:-4]
            
            # if file != 'GM_D18_D25_D32_D37.csv':
            PD_predictions = pd.read_csv(f'{in_dir}/{days_set}/{perc_case}/{file}', index_col=0)
            # else:
            #     PD_predictions = pd.read_csv(f'{wd}/input/Predictions/D18_D25_D32_D37/GM_D18_D25_D32_D37.csv', index_col=0)

            PD_predictions = PD_predictions.index.to_list()            
    
            knownDTIs_l = {}
            drugable = 0
            # not_drugable = 0
            for gene in PD_predictions:
                if gene2drug_ds[gene] != []:
                    knownDTIs_l[gene] = gene2drug_ds[gene]
                    drugable+=1
                else:
                    knownDTIs_l[gene] = []
                    # not_drugable+=1
            try:
                perc_drugable = drugable/len(PD_predictions)*100
                print(drugable)
                print(perc_drugable)
                # print(drugable)
            except:
                print('No genes')
             
            # if file != 'GM_D18_D25_D32_D37.csv':    
            sd = f'{wd}/output/{days_set}/{perc_case}/{save_file}'
            # else:
            #     sd = f'output/D18_D25_D32_D37/{save_file}'
            if not os.path.exists(sd):
                os.makedirs(sd)            
          
            drugs_unique = []
            with open(f'{sd}/{save_file}_druggability.txt', 'w') as f:
                for gene in knownDTIs_l:
                    f.write(f'{gene}\t')
                    # try:
                    n = len(knownDTIs_l[gene])-1
                    if n != -1:
                        # print(n)
                        for i, drug in enumerate(knownDTIs_l[gene]):
                            # print(knownDTIs_l[gene], i,n)
                            drugs_unique.append(drug)
                            name = drug2name[drug]#list(set(DTI[DTI['DrugBank ID'] == drug]['DrugBank Name']))[0]
                            if i != n:
                                f.write(f'{name}, ')
                            else:
                                f.write(f'{name}\n')
                        # except:
                    else:
                        f.write('\n')
          
            count = Counter(drugs_unique)
            count = count.most_common()
            count = [list(x) for x in count]
            # print(count)
            
            count_drug2gene = {}
            for i in range(len(count)):
                pair = count[i]
                drug = pair[0]
                for gene in knownDTIs_l:
                    if drug in knownDTIs_l[gene]:
                        # print(drug, gene)
                        if drug2name[drug] not in count_drug2gene:
                            count_drug2gene[drug2name[drug]] = [gene]
                        else:
                            count_drug2gene[drug2name[drug]].append(gene)
                                  
            with open(f'{sd}/{save_file}_drugs.txt', 'w') as f:
                for i in range(len(count)):
                    pair = count[i]
                    name = drug2name[pair[0]]
                    f.write(f'{name}\t{pair[1]}\t{count_drug2gene[name]}\n')           
    
                    
            # Drugs_dict[perc_case][save_file] = count[:10]
    
    
            bh_drugs, mne_drugs, eg_drugsannot, _, enriched_genes_indrugs = go_enrichment([PD_predictions], gene2drug_ds, drug_count)
            enr_drugs_id =  bh_drugs[0]
            enr_drugs_id = Sort(enr_drugs_id)
            # print(enr_drugs_id)
            enr_drugs_id = [x[0] for x in enr_drugs_id]
            
            enr_drugs = bh_drugs[0]
            enr_drugs = [[drug2name[x[0]],x[1]] for x in enr_drugs]
            enr_drugs = Sort(enr_drugs)
            enriched_genes_indrugs = enriched_genes_indrugs[0]
            # print(enr_drugs, eg_drugsannot, enriched_genes_indrugs)   
            
            # enriched drugs that target at least 2 genes
            enr_drugs_2genes = [drug2name[x] for x in enr_drugs_id if len(count_drug2gene[drug2name[x]])>1]
            top6drugs = ['Artenimol', 'Geldanamycin', 'Rifabutin', 'SNX-5422', 'Tanespimycin', 'Phenethyl Isothiocyanate', 'Copper']
            
            with open(f'{sd}/{save_file}_EnrDrugs.txt', 'w') as f:
                for pair in enr_drugs:
                    f.write(f'{pair[0]}\t{pair[1]}\n')      
            with open(f'{sd}/{save_file}_EnrGenes.txt', 'w') as f:
                for gene in enriched_genes_indrugs:
                    g_drugs = gene2drug_ds[gene]
                    g_drugs = [drug2name[drug] for drug in g_drugs]
                    f.write(f'{gene}\t{g_drugs}\n')     

            
            # Interesting pathways
            for file_Rpgenes in os.listdir(f'{wd}/input/RPBPSubgroupsGenes/{days_set}/{perc_case}'):
                print(file_Rpgenes)
                RPgene2Drug = {}
                used_genes = set()
                with open(f'{wd}/input/RPBPSubgroupsGenes/{days_set}/{perc_case}/{file_Rpgenes}', 'rb') as f:
                    loaded_dict = pickle.load(f)
                for root in loaded_dict:
                    root_dict = {}
                    for subgroup in loaded_dict[root]:
                        subgroup_dict = {}
                        for gene in loaded_dict[root][subgroup]:
                            try:
                                drugs_gene = gene2drug[gene]
                                enrdrugs_gene = list(set(drugs_gene)&set(enr_drugs_id))
                                if len(enrdrugs_gene) != 0:
                                    subgroup_dict[gene]=[drug2name[x] for x in enrdrugs_gene if drug2name[x] in top6drugs]
                                    used_genes.add(gene)
                            except Exception:
                                pass
                        root_dict[subgroup] = subgroup_dict
                    RPgene2Drug[root] = root_dict
                # print(RPgene2Drug)
                # print(used_genes)
                
                with open(f'{sd}/{file_Rpgenes[:-4]}_Path2Gene2Drug.txt', 'w') as f:
                    for rootpath in RPgene2Drug:
                        f.write(f'{rootpath}\n')
                        for path in RPgene2Drug[rootpath]:
                            f.write(f'{path}\n')
                            for gene in RPgene2Drug[rootpath][path]:
                                f.write(f'{gene}\t{RPgene2Drug[rootpath][path][gene]}\n')
                        f.write('\n')
                                

                                    
                            
                        
                    



                    
# print('DB11093', drug_count['DB11093']) #Calcium citrate
# print('DB11348', drug_count['DB11348']) #Calcium Phosphate

# names = {'R-HSA-72613':'Eukaryotic Translation Initiation','R-HSA-392499':'Metabolism of proteins'
#          ,'GO:0006364':'rRNA processing','GO:0017148':'negative regulation of translation', 'GO:0031640':'killing of cells of another organism'}

# with open(f'{wd}/input/RP_representative_terms_genes.pkl', 'rb') as f:
#     loaded_dict = pickle.load(f)

# RepGenes_drugs_RP = {names[key]:{} for key in loaded_dict}
# for path in loaded_dict:
#     for gene in loaded_dict[path]:
#         try:
#             drugs = gene2drug[gene]
#             drugs_names = [drug2name[drug] for drug in drugs]
#             RepGenes_drugs_RP[names[path]][gene]=drugs_names
#         except Exception:
#             pass

# genes_paths = [list(RepGenes_drugs_RP[path].keys()) for path in RepGenes_drugs_RP]
# print(set.intersection(*map(set,genes_paths)))
# print(set(genes_paths[1])-set(genes_paths[0]))

# for path in RepGenes_drugs_RP:
#     print(path,RepGenes_drugs_RP[path].keys())


# with open(f'{wd}/input/BP_representative_terms_genes.pkl', 'rb') as f:
#     loaded_dict = pickle.load(f)

# RepGenes_drugs_BP = {names[key]:{} for key in loaded_dict}
# for path in loaded_dict:
#     for gene in loaded_dict[path]:
#         try:
#             drugs = gene2drug[gene]
#             drugs_names = [drug2name[drug] for drug in drugs]
#             RepGenes_drugs_BP[names[path]][gene]=drugs_names
#         except Exception:
#             pass

# for path in RepGenes_drugs_BP:
#     print(path,RepGenes_drugs_BP[path].keys())


