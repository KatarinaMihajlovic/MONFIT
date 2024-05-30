# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:40:06 2024

@author: Katarina
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sd = 'output/GMM_RMSE'     
cell_conds = ['ND_8', 'ND_18', 'ND_25', 'ND_32', 'ND_37', 'WT_8', 'WT_18', 'WT_25', 'WT_32', 'WT_37']

ccs_RI = []           # run1_RS.append(numClst2_RS)
for cc in cell_conds:
    print(cc)
    cc_RS_df = pd.read_csv(f'{sd}/{cc}_GMM_RMSE.csv',index_col=0, header=0)
    labels = list(cc_RS_df.index)
    labels = [str(x) for x in labels]
    cc_RS_df = cc_RS_df.rename(columns=dict(zip(cc_RS_df.columns, labels)))
    
    # Get the color palette
    palette = sns.color_palette('pastel', 5)        
    # Convert the palette to a list
    colors = palette.as_hex()    
    lut = dict(zip(set(labels), colors))
    row_colors = pd.DataFrame(labels)[0].map(lut)

    
    # fig, ax = plt.subplots(figsize=(8, 8))
    sns.clustermap(cc_RS_df, cmap='viridis', col_cluster=False, row_cluster=False, row_colors=row_colors, col_colors=row_colors)
    # handles = [Patch(facecolor=lut[name]) for name in lut]
    # plt.legend(handles, lut, title='Number of Clusters', fontsize = 18, title_fontsize=18,
    #            bbox_to_anchor=(0.2, 0.4), bbox_transform=plt.gcf().transFigure, loc='upper right')
    # plt.tight_layout()  
    plt.show()
    
    # median, stdv of RS
    upper_triangle_values = np.triu(cc_RS_df.values, k=1)
    # Convert upper triangle values to a 1D array
    upper_triangle_values = upper_triangle_values[np.triu_indices(cc_RS_df.values.shape[0])]
    upper_triangle_values = upper_triangle_values[upper_triangle_values != 0]
    ccs_RI.append(upper_triangle_values)
    
    mean = np.mean(upper_triangle_values)
    stdv = np.std(upper_triangle_values)
    print(mean, stdv)

print('stats for RMSE across all CCs')        
ccs_RI =  [item for row in ccs_RI for item in row]

mean = np.mean(ccs_RI)
stdv = np.std(ccs_RI)
print(mean, stdv)