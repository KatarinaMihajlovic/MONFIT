# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:05:03 2024

@author: Katarina
"""

# IMPORTANT NOTE
# This script requires precomputed NMTFs with different weighting factors.

import pandas as pd
import os
from sklearn.metrics.pairwise import euclidean_distances 
import numpy as np
from sklearn.metrics import mean_squared_error
import math

def write_txt(file_path, name, list_l):
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    with open(f'{file_path}/{name}', 'w') as f: 
        for x in list_l:
            f.write(f'{x}\n')




work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

G1S5_dir = 'step2_NMTF_WeightedNets/output'
avg_score_df = pd.read_csv(f'{work_dir}/input/Enrichment_Rank.csv', header = 0, index_col=0)
best_combs = list(avg_score_df.index)[:100]



sd = f'{work_dir}/output/GMM_RMSE'
try:
    os.makedirs(sd)
except FileExistsError:
    # directory already exists
    pass

cell_conds = ['ND_8', 'ND_18', 'ND_25', 'ND_32', 'ND_37', 'WT_8', 'WT_18', 'WT_25', 'WT_32', 'WT_37']
for cc in cell_conds:
    cc_GMM_RSE = np.empty((len(best_combs), len(best_combs)))   
    
    for k1k2 in os.listdir(f'{G1S5_dir}/{cc}'):
        for i in range(len(best_combs)):
            comb1 = best_combs[i]
            G1 = pd.read_csv(f'{G1S5_dir}/{cc}/{k1k2}/{comb1}/G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
            S5 = pd.read_csv(f'{G1S5_dir}/{cc}/{k1k2}/{comb1}/S_EXP.csv', header=0, index_col=0, delimiter='\t')
    
            G1_vals = G1.to_numpy()
            S5_vals = S5.to_numpy()
            G1_S5 = G1_vals.dot(S5_vals)
           
            # Normalized Euclidean distance       
            EuclidDist_all = euclidean_distances(G1_S5)
            norm = np.linalg.norm(EuclidDist_all) #, axis = 0)
            EuclidDist_all_normalized_1 = EuclidDist_all/norm
            for j in range(i, len(best_combs)):
                comb2 = best_combs[j]
                G1 = pd.read_csv(f'{G1S5_dir}/{cc}/{k1k2}/{comb2}/G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
                S5 = pd.read_csv(f'{G1S5_dir}/{cc}/{k1k2}/{comb2}/S_EXP.csv', header=0, index_col=0, delimiter='\t')
    
                G1_vals = G1.to_numpy()
                S5_vals = S5.to_numpy()
                G1_S5 = G1_vals.dot(S5_vals)
               
                # Normalized Euclidean distance       
                EuclidDist_all = euclidean_distances(G1_S5)
                norm = np.linalg.norm(EuclidDist_all) #, axis = 0)
                EuclidDist_all_normalized_2 = EuclidDist_all/norm                
                    
                mse = mean_squared_error(EuclidDist_all_normalized_1, EuclidDist_all_normalized_2)
                rmse = math.sqrt(mse)  
                cc_GMM_RSE[i][j] = rmse
                cc_GMM_RSE[j][i] = rmse
    cc_GMM_RSE_df = pd.DataFrame(cc_GMM_RSE,columns=best_combs,index=best_combs)
    cc_GMM_RSE_df.to_csv(f'{sd}/{cc}_GMM_RMSE.csv')                
                    
os.chdir(work_dir) 
                                
   