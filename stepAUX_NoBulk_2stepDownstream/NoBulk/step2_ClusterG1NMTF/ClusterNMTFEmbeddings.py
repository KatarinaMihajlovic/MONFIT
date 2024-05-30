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


work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

cell_conds = ['ND_8', 'ND_18', 'ND_25', 'ND_32', 'ND_37', 'WT_8', 'WT_18', 'WT_25', 'WT_32', 'WT_37']
for cc in cell_conds:   
    k1k2 = os.listdir(f'step1_NMTF/output/{cc}')[0]
    for comb in os.listdir(f'step1_NMTF/output/{cc}/{k1k2}'):
        save_dir = f'{work_dir}/output/{cc}/{comb}'
        try:
            os.makedirs(save_dir)
        except FileExistsError:
            # directory already exists
            pass

        save_file = f'G1_clsts_{comb}'
        if len(os.listdir(save_dir)) != 10:
            try:
                G1_df = pd.read_csv(f'step1_NMTF/output/{cc}/{k1k2}/{comb}/G1_with_headers.csv', header=0, index_col=0, sep='\t')
                G1_NumClusters = len(list(G1_df.columns))
                genes = list(G1_df.index)
                G1 = G1_df.values
        
                for run in range(km_runs):
                    clustersG1_kmean = kMeans(G1, genes, G1_NumClusters, km_runs)
                    with open(f'{save_dir}/{run}_{save_file}', 'wb') as handle:
                        pickle.dump(clustersG1_kmean, handle)
            except Exception as e:
                print(e)


os.chdir(work_dir)     
