# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 20:00:30 2021

@author: kmihajlo
"""

import pickle, os


def Sort(sub_li):
  
    sub_li.sort(key = lambda x: x[1], reverse=True)
    return sub_li



### Main Code

with open('input/PDgenes_DGN_ALL.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)  


for cc in os.listdir('output'):
    for file in os.listdir(f'output/{cc}'):
        if 'indices' not in file and 'pkl' in file and 'InitPreds' not in file:
            print(cc,file)
            Init_preds = []
            
            print(f'output/{cc}/{file}')
            with open(f'output/{cc}/{file}', 'rb') as handle:
                ClustsPDEnrich = pickle.load(handle)
            
            i = 0
            for clst in range(len(ClustsPDEnrich)):
                PD_gene_Predictions = []
                EnrClust = ClustsPDEnrich[clst]
                Litvalid_genes = []
                for gene in EnrClust:
                    if gene not in PD_genes:
                        PD_gene_Predictions.append(gene)
                        i+=1
                Init_preds.append(PD_gene_Predictions)
            Init_preds = [j for i in Init_preds for j in i]
            print(len(Init_preds))
                
            with open(f'output/{cc}/{cc}_InitPreds.txt', 'w') as f:
                for gene in Init_preds:
                    f.write(f'{gene}\n')
            with open(f'output/{cc}/{cc}_InitPreds.pkl', 'wb') as handle:
                pickle.dump(Init_preds, handle)


                        
# with open('input/Clusters/Control_D21/kMeans_PDEnrichClusts.pkl', 'rb') as handle:
#     dsa = pickle.load(handle)  
                        
                        