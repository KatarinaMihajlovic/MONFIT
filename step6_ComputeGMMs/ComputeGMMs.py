# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 12:38:35 2021

@author: kmihajlo
"""

import pandas as pd
import os, sys
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np

wd = str(sys.argv[1])


for cc in os.listdir(f'{wd}/input/NMTF_G1s'): #C1 = cell condition 1
    print(cc)
    G1 = pd.read_csv(f'{wd}/input/NMTF_G1s/{cc}/{cc}_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
    S5 = pd.read_csv(f'{wd}/input/NMTF_G1s/{cc}/{cc}_S_EXP.csv', header=0, index_col=0, delimiter='\t')

    G1_vals = G1.to_numpy()
    S5_vals = S5.to_numpy()
    G1_S5 = G1_vals.dot(S5_vals)
    
    # Normalized Euclidean distance       
    EuclidDist_all = euclidean_distances(G1_S5)
    norm = np.linalg.norm(EuclidDist_all) #, axis = 0)
    EuclidDist_all_normalized = EuclidDist_all/norm

    sd = f'{wd}/output'
    if not os.path.exists(sd):
        os.makedirs(sd) 
    save_file = f'{sd}/{cc}_PWEucDist'
    np.save(save_file, EuclidDist_all_normalized)  
    
   