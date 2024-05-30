# -*- coding: UTF-8 -*-

''' Non-negative matrix tri-factorization (numpy)'''
# Author: kmihajlo
import os
import networkx as nx
import numpy as np
import Network_Matrices as nm
import Matrix_Factorization as mf
import time
import sys
import pandas as pd
import matplotlib.pyplot as plt

os.environ['OMP_NUM_THREADS'] = '8' # export OMP_NUM_THREADS=4
os.environ['OPENBLAS_NUM_THREADS'] = '8' # export OPENBLAS_NUM_THREADS=4 
os.environ['MKL_NUM_THREADS'] = '8' # export MKL_NUM_THREADS=6
os.environ['VECLIB_MAXIMUM_THREADS'] = '8' # export VECLIB_MAXIMUM_THREADS=4
os.environ['NUMEXPR_NUM_THREADS'] = '8' # export NUMEXPR_NUM_THREADS=6




print('\n\x1b[1;37;44m ############################################# \x1b[0m')
print('\x1b[1;37;44m #                                           # \x1b[0m')
print('\x1b[1;37;44m # ISF Framework                             # \x1b[0m')
print('\x1b[1;37;44m #                                           # \x1b[0m')
print('\x1b[1;37;44m ############################################# \x1b[0m')

epsilon = 1e-15

cell_cond = str(sys.argv[1])
k1 = int(sys.argv[2])
k2 = int(sys.argv[3])


E_w = 10.0
PPI_w = 1.0
COEX_w = 1.0
GI_w = 0.0
MI_w = 10.0

print(cell_cond)

# P_1.0_C_10.0_G_0.1_M_0.0_E_0.1
# PPI_ws = [0, 0.1, 1, 10]
# COEX_ws = [0, 0.1, 1, 10]
# GI_ws = [0, 0.1, 1, 10]
# MI_ws = [0, 0.1, 1, 10]
# E_ws = [0.1, 1, 10]
   
# ND_18 0.1 1 10 0.1 0
# cell_cond = 'WT_18'
# E_w = 0.1
# PPI_w = 10.0
# COEX_w = 0.1
# GI_w = 1.0
# MI_w = 1.0

# WT_18/P_10.0_C_0.1_G_0.0_M_1.0_E_0.1


# for cell_cond in os.listdir('input'):  
#     print(cell_cond)    
print('---- Loading Expression matrix')
EXP = pd.read_csv(f'input/{cell_cond}/E_{cell_cond}_normalized.csv', index_col=0)
g,c = EXP.shape
# k1 = int(np.sqrt(g/2))
# k2 = int(np.sqrt(c/2))   
k1k2_folder = f'output/{cell_cond}/k1_{k1}_k2_{k2}'     
save_dir = f'{k1k2_folder}/P_{PPI_w}_C_{COEX_w}_G_{GI_w}_M_{MI_w}_E_{E_w}' 
print(f'P_{PPI_w}_C_{COEX_w}_G_{GI_w}_M_{MI_w}_E_{E_w}')

try:
    os.makedirs(save_dir)
except FileExistsError:
    # directory already exists
    pass
    
print(len(os.listdir(save_dir)))
if len(os.listdir(save_dir)) == 0:

    nodes = list(EXP.index)
    node2ind = {}
    for n in range(len(nodes)):
        node2ind[nodes[n]]=n
    # EXP = EXP.loc[nodes]
    EXP = EXP.values
    EXP = E_w*EXP
    EXP+=epsilon
              
    print('---- Loading PPI network')
    # PPI, nodes, node2ind = nm.Load_Network_weighted(f'input/{cell_cond}/PPI_ppTW_{cell_cond}.edgelist') # nodes = gene nodes, node2ind = gene node to index
    if PPI_w != 0:
        PPI = nx.read_weighted_edgelist(f'input/{cell_cond}/PPI_ppTW_{cell_cond}.edgelist')
        PPI_mat = nm.Make_Adj_weighted(PPI, nodes, node2ind)
        PPI_mat = PPI_w*PPI_mat
    
    if COEX_w != 0:
        print('---- Loading COEX network')
        COEX = nx.read_edgelist(f'input/{cell_cond}/COEX_{cell_cond}.edgelist')
        COEX_mat = nm.Make_Adj(COEX, nodes, node2ind)
        COEX_mat = COEX_w*COEX_mat
    
    if GI_w != 0:
        print('---- Loading GI network')
        GI = nx.read_edgelist(f'input/{cell_cond}/GI_{cell_cond}.edgelist')
        GI_mat = nm.Make_Adj(GI, nodes, node2ind)
        GI_mat = GI_w*GI_mat
    if MI_w != 0:
        print('---- Loading MI network')
        MI = nx.read_weighted_edgelist(f'input/{cell_cond}/MI_MetAbundWeight_{cell_cond}.edgelist')
        MI_mat = nm.Make_Adj_weighted(MI, nodes, node2ind)
        MI_mat = MI_w*MI_mat
            
    
    
    with open(f'output/Geneslist_{cell_cond}.csv', 'w') as f:
        f.write('genes' + '\n')
        for gene in nodes:
            f.write(f'{gene}\n')
                    
            
    #start NMTF    
    start = time.time()
    print(k1, k2)
    
    featuresG1 = [str(i) for i in range(k1)]
    featuresG2 = [str(i) for i in range(k2)]
    
    
    
    
    if len(os.listdir(save_dir)) == 0:               
        MOLs = []
        used_nets = []
        if PPI_w!=0:
            MOLs.append(PPI_mat)
            used_nets.append('PPI')
        if COEX_w!=0:
            MOLs.append(COEX_mat)
            used_nets.append('COEX')
        if GI_w!=0:
            MOLs.append(GI_mat)
            used_nets.append('GI')
        if MI_w!=0:
            MOLs.append(MI_mat)
            used_nets.append('MI')
        if PPI_w == COEX_w == GI_w == MI_w == 0:
            used_nets = ['noMOLs']
            MOLs = [np.zeros((g,g))+epsilon]
        print(used_nets)
            
        if PPI_w == 0 and COEX_w == 0 and GI_w == 0 and MI_w == 0:
            Solver = mf.NMTF_basic(max_iter=1000, verbose = 10)
            G1, G2, S_EXP, OBJ_fig = Solver.Solve_MUR(EXP, k1, k2, init='SVD')
        else:
            Solver = mf.PD_SSNMTF(max_iter=1000, verbose = 10)
            G1, G2, S_Mol, S_EXP, OBJ_fig = Solver.Solve_MUR(MOLs, EXP, k1, k2, init='SVD')
            for i in range(len(used_nets)):
                nm.Save_Matrix_Factor(S_Mol[i], save_dir  + '/' + 'S_' + used_nets[i] + '.csv', featuresG1, featuresG1)	
      
            
        OBJ_fig.tight_layout()    
        OBJ_fig.savefig(save_dir + '/' + 'OBJ', dpi = 350)
        plt.close(OBJ_fig)
        
        print('saving matrices')
        nm.Save_Matrix_Factor(G1, save_dir + '/' +'G1_with_headers.csv', nodes, featuresG1)
        nm.Save_Matrix_Factor(G2, save_dir + '/'  +'G2_with_headers.csv', nodes, featuresG2)
        nm.Save_Matrix_Factor(S_EXP, save_dir  + '/' + 'S_EXP.csv', featuresG1, featuresG2)	
    
    end = time.time()
    run_time = (end - start)/60
    print('Runtime of the program is ' + str(run_time))
    
       

    
       

