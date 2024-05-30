# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 12:38:35 2021

@author: kmihajlo
"""

import pandas as pd
import os, pickle
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances, manhattan_distances
import numpy as np
import matplotlib.pyplot as plt


def write_txt(file_path, name, list_l):
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    with open(f'{file_path}/{name}', 'w') as f: 
        for x in list_l:
            f.write(f'{x}\n')


ccs = os.listdir(f'input/NMTF_matrices')

for cci in range(len(ccs)): #C1 = cell condition 1
    cc = ccs[cci]
    G1 = pd.read_csv(f'input/NMTF_matrices/{cc}/{cc}_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
    S5 = pd.read_csv(f'input/NMTF_matrices/{cc}/{cc}_S_EXP.csv', header=0, index_col=0, delimiter='\t')

    G1_vals = G1.to_numpy()
    S5_vals = S5.to_numpy()
    G1_S5 = G1_vals.dot(S5_vals)


    # Normalized Euclidean distance
    EuclidDist_all = euclidean_distances(G1_S5)
    feature_norms = np.linalg.norm(EuclidDist_all) #, axis = 0)
    EuclidDist_all_normalized = EuclidDist_all/feature_norms

    
    file_path = 'output/EuclideanDistance'
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    save_file = f'{file_path}/{cc}_PWEucDist'
    np.save(save_file, EuclidDist_all_normalized)  
    

    with open(f'output/{cc}_Genes.pkl', 'wb') as handel:
        pickle.dump(G1.index, handel)    


