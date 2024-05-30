# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 16:49:02 2024

@author: Katarina
"""

import pickle, os, sys
from sklearn.cluster import KMeans
import pandas as pd


def kMeans(M, nodes, n_clusters, km_runs):
    M = M.tolist()
    nodes2coordinates = dict(zip(nodes, M))  
    Cluster_belonging = {k: [] for k in range(n_clusters)}
    
    try:
        kMeans_model = KMeans(n_clusters=n_clusters, init = 'random', n_init=km_runs).fit(M)
    except Exception:
        print(M)
    KMeans_labels = list(kMeans_model.labels_)
    
    for cluster_index in range(len(KMeans_labels)):
        cluster = KMeans_labels[cluster_index]
        node_coords = M[cluster_index]
        for node, coordinates in nodes2coordinates.items():
            if node_coords == coordinates:
                Cluster_belonging[cluster].append(node)
    
    #print(Cluster_belonging)
    Cluster_belonging_list = []
    for _, values in Cluster_belonging.items():
        Cluster_belonging_list.append(values)
    #print(Cluster_belonging_list)

    return Cluster_belonging_list

km_runs = 10


E_w = float(sys.argv[1])
PPI_w = float(sys.argv[2])
COEX_w = float(sys.argv[3])
GI_w = float(sys.argv[4])
MI_w = float(sys.argv[5])
wd = str(sys.argv[6])

# cell_cond = 'ND_18'
# E_w = 0.1
# PPI_w = 1.0
# COEX_w = 10.0
# GI_w = 0.1
# MI_w = 10.0


in_dir = f'{wd}/input/NMTF_G1s' #'/P_{PPI_w}_C_{COEX_w}_G_{GI_w}_M_{MI_w}_E_{E_w}' 
ccs = os.listdir(in_dir)
comb = f'P_{PPI_w}_C_{COEX_w}_G_{GI_w}_M_{MI_w}_E_{E_w}'

for cc in ccs:
    save_dir = f'{wd}/output/{cc}/{comb}' 
    save_file = f'G1_clsts_{comb}.pkl'       
    
    try:
        os.makedirs(save_dir)
    except FileExistsError:
        # directory already exists
        pass
    
    if len(os.listdir(save_dir)) != 10:
        try:
            G1_df = pd.read_csv(f'{in_dir}/{cc}/{comb}/{cc}_G1_with_headers.csv', header=0, index_col=0, sep='\t')
            G1_NumClusters = len(list(G1_df.columns))
            genes = list(G1_df.index)
            G1 = G1_df.values
            print(cc, f'P_{PPI_w}_C_{COEX_w}_G_{GI_w}_M_{MI_w}_E_{E_w}')
    
            for run in range(km_runs):
                clustersG1_kmean = kMeans(G1, genes, G1_NumClusters, km_runs)
                with open(f'{save_dir}/{run}_{save_file}', 'wb') as handle:
                    pickle.dump(clustersG1_kmean, handle)
        except Exception as e:
            print(e)
    
