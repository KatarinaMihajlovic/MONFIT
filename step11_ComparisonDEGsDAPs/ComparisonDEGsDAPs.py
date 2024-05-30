# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:37:32 2023

@author: kmihajlo
"""

import pandas as pd
import os, sys
import matplotlib.pyplot as plt
import venn
import warnings
warnings.filterwarnings("ignore")

def write_txt(list_l, file_path, name):
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    with open(f'{file_path}/{name}', 'w') as f: 
        for x in list_l:
            f.write(f'{x}\n')    
 

def fc_preds(PD_predictions, DEGs_ds):
    DEG_FC = []
    for pred in PD_predictions:
        if pred in list(DEGs_ds.index):
            DEG_FC.append(DEGs_ds.loc[pred]['absFCpower'])
        else:
            DEG_FC.append(0)
    return DEG_FC

        
'''MAIN'''
wd = str(sys.argv[1])

    
DEGs = pd.read_excel(f'{wd}/input/DEG_WT_vs_ND_0to57_SUMMARY.xlsx')
DEGs = DEGs.drop(['0','57'], axis=1)
DEGs = DEGs.rename(index=dict(zip(DEGs.index, DEGs['Gene']))).drop(['Gene','absFCpower'],axis=1)
DEGs = DEGs.loc[(abs(DEGs).sum(axis=1) != 0)]


DAPs = pd.read_excel(f'{wd}/input/DAP_WT_VS_ND_0to57_SUMMARIZED.xlsx')
DAPs = DAPs.drop(['0','57'], axis=1)
DAPs = DAPs.rename(index=dict(zip(DAPs.index, DAPs['Gene']))).drop(['Gene','absFCpower'],axis=1)
DAPs = DAPs.loc[(abs(DAPs).sum(axis=1) != 0)]


in_dir = f'{wd}/input/Predictions'

for days_set in os.listdir(in_dir): 
    days = days_set.split('_')
    days = [x[1:] for x in days]
    sd = f'{wd}/output/{days_set}'
    if not os.path.exists(sd):
        os.makedirs(sd) 

    # DEGs
    DEGs_ds = DEGs[days].copy()
    DEGs_ds = DEGs_ds.loc[(abs(DEGs_ds).sum(axis=1) != 0)]     
    thresh_DEG = 0.5#elbowp #-float('inf') #0.3 
    DEGs_l = []   
    for day in days:
        genes = list(DEGs_ds[abs(DEGs_ds[day]) > thresh_DEG].index)
        DEGs_l.append(genes)
    DEGs_l = list(set().union(*DEGs_l))
                    
    DEGs_ds = DEGs_ds.loc[DEGs_l]
    DEGs_ds['absFCpower'] = abs(DEGs_ds).mean(axis=1)
    DEGs_ds = abs(DEGs_ds).sort_values(by=['absFCpower', *days], ascending=False)
    DEGs_ds.to_csv(f'{sd}/DEGs_{days_set}_{thresh_DEG}FC.csv')

    # DAPs              
    DAPs_ds = DAPs[days].copy()
    DAPs_ds = DAPs_ds.loc[(abs(DAPs_ds).sum(axis=1) != 0)]
    thresh_DAPs = 1#elbowp #-float('inf') #0.3 
    DAPs_l = []   
    for day in days:
        genes = list(DAPs_ds[abs(DAPs_ds[day]) > thresh_DAPs].index)
        DAPs_l.append(genes)
    DAPs_l = list(set().union(*DAPs_l))
        
    DAPs_ds = DAPs_ds.loc[DAPs_l]    
    DAPs_ds['absFCpower'] = abs(DAPs_ds).mean(axis=1)
    DAPs_ds = abs(DAPs_ds).sort_values(by=['absFCpower', *days], ascending=False)
    DAPs_ds.to_csv(f'{sd}/DAPs_{days_set}_{thresh_DAPs}FC.csv')
        
    print(len(DEGs_l), len(DAPs_l))
    # days_set = 'D8_D18_D25_D32_D37'
    for perc_case in os.listdir(f'{in_dir}/{days_set}'):
        # perc_case = '1.157perc'
        for file in os.listdir(f'{in_dir}/{days_set}/{perc_case}'):
            # file = 'GM_Length_D8_D18_D25_D32_D37_1.157_PDpreds.csv'
            save_file = file[10:-4]
            print(save_file)
            sd = f'{wd}/output/{days_set}/{perc_case}/DEG{thresh_DEG}_DAP{thresh_DAPs}'

            if not os.path.exists(sd):
                os.makedirs(sd) 
                
                
            PD_predictions_df = pd.read_csv(f'{in_dir}/{days_set}/{perc_case}/{file}', index_col=0)
            PD_predictions = PD_predictions_df.index.to_list()    
            
            # ribosomal = []
            # for gene in PD_predictions:
            #     if 'RPL' == gene[:3] or 'RPS' == gene[:3]:
            #         ribosomal.append(gene)
            # print(len(ribosomal))
            
            
            PDpreds_DEGs = list(set(DEGs_l) & set(PD_predictions))
            PDpreds_DAPs = list(set(DAPs_l) & set(PD_predictions))

            # Venn
            genes = [set(PD_predictions), set(DEGs_l),set(DAPs_l)]
            gene_sets = ['MONFIT', 'DEGs','DAPs']
            labels = venn.generate_petal_labels(genes)#, fill=['number', 'logic'])
            fig, ax = venn.venn3(labels, names=gene_sets, fontsize =24, figsize= (9, 12))
            ax.set_title('Genes', fontsize = 28, y=0.9) # Gene Overlap between\nPD predictions, DEGs and DAPs'   
            ax.get_legend().remove()
            # leg = ax.legend(gene_sets, loc='lower left', bbox_to_anchor=(0, -0.14), fancybox=True, fontsize = 22)
            plt.savefig(f'{sd}/{save_file}_GeneOverlap_DEG{thresh_DEG}_DAP{thresh_DAPs}.jpg',  dpi = 350, format='jpg', bbox_inches='tight')	
            plt.show()


            # saving gene groups
            sd_gg = f'{sd}/GeneGroups'
            if not os.path.exists(sd_gg):
                os.makedirs(sd_gg) 
                
                
            # geneU = set().union(*genes)
            # set.intersection(*map(set,d))
            Preds_only = list(set(PD_predictions)-set(DAPs_l)-set(DEGs_l))
            DEGs_only = list(set(DEGs_l)-set(PD_predictions)) #-set(DAPs_l)
            DAPs_only = list(set(DAPs_l)-set(PD_predictions)) #-set(DEGs_l)
            Intersection_all = list(set(PD_predictions)&set(DAPs_l)&set(DEGs_l))

            write_txt(Preds_only, sd_gg, 'MONFIT_only.txt')
            write_txt(DEGs_only, sd_gg, 'DEGs_only.txt')
            write_txt(DAPs_only, sd_gg, 'DAPs_only.txt')
            write_txt(Intersection_all, sd_gg, 'MultiValidated.txt')

            
            PD_predictions_DEG_FC = fc_preds(PD_predictions, DEGs_ds)
            PD_predictions_DAP_FC = fc_preds(PD_predictions, DAPs_ds)          
            PD_preds_FC = pd.DataFrame([PD_predictions,PD_predictions_DEG_FC,PD_predictions_DAP_FC], index=['Preds','DEGs','DAPs']).T
                


            # top 30 
            top30Preds = PD_predictions[:30]
            Pred_DEG = list(set(top30Preds)-set(DAPs_l)&set(DEGs_l))
            Pred_DAP = list(set(top30Preds)&set(DAPs_l)-set(DEGs_l))
            Pred_DEG_DAP = list(set(top30Preds)&set(DAPs_l)&set(DEGs_l))
            print(Pred_DEG, Pred_DAP, Pred_DEG_DAP)
          
            
            
            
            
            
    