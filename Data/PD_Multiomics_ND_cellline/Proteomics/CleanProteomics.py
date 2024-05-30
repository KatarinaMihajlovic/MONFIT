# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:30:53 2023

@author: kmihajlo
"""

import pyreadr
import pandas as pd
import numpy as np 

def avg_replicates(df, objects, cc):
    abundances = []
    for obj in objects:
        abundance_reps = np.array(df[obj])
        abundance_reps = np.append(abundance_reps,[0,0])
        zeros = len(abundance_reps) - np.count_nonzero(abundance_reps)
        
        if zeros == 0:
            abundance = np.mean(abundance_reps)
        elif zeros <= len(abundance_reps)/2:
            abundance_reps_no_0 = abundance_reps[abundance_reps!=0]
            abundance = np.mean(abundance_reps_no_0)
        else:
            abundance = 0
        abundances.append(abundance)
    
    df_avg = pd.DataFrame(abundances,columns=[f'{cc}'],index=objects)
    return df_avg
        


### PROTEOMICS
# avreage over replicates, remove duplicate proteins a nd log scale the abundances
result = pyreadr.read_r('input/P_spread_WT_vs_ND_8to37.RData')
Proteomics = result[None]
result = pyreadr.read_r('input/ProteinToGene.RData')
protein2gene = result[None]

protein2gene_dict = dict(zip(protein2gene.index.to_list(), protein2gene.Gene))
Proteomics = Proteomics.rename(columns=protein2gene_dict)
gene_proteins = Proteomics.columns.to_list()[3:]

days = set(Proteomics['Day'])
days = [int(x) for x in days]
types_cc = set(Proteomics['Condition'])

ccs = []
Proteomics_avg_log = []

# avg across replicates and split into conditinos 
for type_cc in types_cc:
    for day in days:
        cc = f'{type_cc}_{day}'
        print(cc)
        ccs.append(cc)
        cc_Prot = Proteomics.loc[(Proteomics['Condition'] == type_cc) & (Proteomics['Day'] == day)]
        genes_p = cc_Prot.columns.to_list()[3:]

        cc_Prot_avg = avg_replicates(cc_Prot, genes_p, cc)   
        cc_Prot_avg['index_column'] = cc_Prot_avg.index
        cc_Prot_avg_nodup = cc_Prot_avg.sort_values(cc, ascending=False).drop_duplicates('index_column').sort_index()
        cc_Prot_avg_nodup = cc_Prot_avg_nodup.drop('index_column', axis =1)
        # cc_Prot_avg_nodup[cc] = cc_Prot_avg_nodup[cc]/cc_Prot_avg_nodup[cc].sum()

        cc_Prot_avg_nodup_log = np.log10(cc_Prot_avg_nodup+1)

        cc_Prot_avg_nodup_vals = cc_Prot_avg_nodup_log[cc].values
        Proteomics_avg_log.append(cc_Prot_avg_nodup_vals)
        genes_p = list(cc_Prot_avg_nodup_log.index)
                


Proteomics_avg_log = np.array(Proteomics_avg_log)
Proteomics_avg_log_df = pd.DataFrame(Proteomics_avg_log, index = ccs, columns = genes_p)  
Proteomics_avg_log_df = Proteomics_avg_log_df.loc[:, Proteomics_avg_log_df.sum(axis=0) > 0]
Proteomics_avg_log_df.to_csv('output/Proteomicslog.csv',index=True)  





