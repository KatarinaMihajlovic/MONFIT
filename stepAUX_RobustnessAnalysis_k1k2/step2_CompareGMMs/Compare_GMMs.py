# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:05:03 2024

@author: Katarina
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 12:38:35 2021

@author: kmihajlo
"""

import pandas as pd
import os, sys
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


cell_conds = ['ND_8', 'ND_18', 'ND_25', 'ND_32', 'ND_37', 'WT_8', 'WT_18', 'WT_25', 'WT_32', 'WT_37']



work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


G1S5_dir = 'step1_NMTF/output'

sd = f'{work_dir}/output/GMM_RMSE'
try:
    os.makedirs(sd)
except FileExistsError:
    # directory already exists
    pass

for cc in cell_conds:
    combs = os.listdir(f'{G1S5_dir}/{cc}')
    best_comb = 'P_1.0_C_1.0_G_0.0_M_10.0_E_10.0'
    cc_GMM_RSE = np.empty((len(combs), len(combs)))   
    
    for i in range(len(combs)):
        G1 = pd.read_csv(f'{G1S5_dir}/{cc}/{combs[i]}/{best_comb}/G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
        S5 = pd.read_csv(f'{G1S5_dir}/{cc}/{combs[i]}/{best_comb}/S_EXP.csv', header=0, index_col=0, delimiter='\t')
    
        G1_vals = G1.to_numpy()
        S5_vals = S5.to_numpy()
        G1_S5 = G1_vals.dot(S5_vals)
       
        # Normalized Euclidean distance       
        EuclidDist_all = euclidean_distances(G1_S5)
        norm = np.linalg.norm(EuclidDist_all) #, axis = 0)
        EuclidDist_all_normalized_1 = EuclidDist_all/norm
        for j in range(i, len(combs)):
            G1 = pd.read_csv(f'{G1S5_dir}/{cc}/{combs[j]}/{best_comb}/G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
            S5 = pd.read_csv(f'{G1S5_dir}/{cc}/{combs[j]}/{best_comb}/S_EXP.csv', header=0, index_col=0, delimiter='\t')
    
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
    cc_GMM_RSE_df = pd.DataFrame(cc_GMM_RSE,columns=combs,index=combs)
    cc_GMM_RSE_df.to_csv(f'{sd}/{cc}_GMM_RMSE.csv')                
                    
os.chdir(work_dir) 
                                
   