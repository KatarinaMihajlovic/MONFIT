# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 14:09:59 2022

@author: kmihajlo
"""
import os, pickle,sys
import numpy as np
import pandas as pd
from scipy.spatial import distance
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kstest
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
from matplotlib.ticker import FormatStrFormatter
from kneed import KneeLocator
 
default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [99, 236, 243, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
] 




def plot_hist_GML(GM_length, v_line, TP, TP_len, title, sd, save_ext):
    fig, ax = plt.subplots(figsize=(11, 9))
    plt.hist(GM_length, bins = 100)#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 
    plt.axvline(x = v_line, color = 'r', label = f'{round(TP,3)}% thresh ({TP_len})')
    ax.set_yscale('log')
    plt.title(title, fontsize=26) 
    ax.set_xlabel('Gene movement length', fontsize=26) #dist_measure_l[dist_meas]
    ax.set_ylabel('frequency', fontsize=26)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=26)
    plt.savefig(f'{sd}/GM_length_{save_ext}_hist.jpg', dpi = 350, format='jpg')  
    plt.show()
    plt.close()      


def PDGs2Background_GM(dist_df, case_compare, TargetGenes, save_dir, stat_meas = 'MWU', PD_Gs = 'PD preds'):
    print(PD_Gs)
    # save_dir = f'output/StageSpecPreds/{sim_measure}/{dist_meas}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)  
        
    GM_genes = []
    GM_BGgenes = []
    for gene in dist_df.index:
        if gene in TargetGenes:    
            GM_genes.append(dist_df.loc[gene,case_compare])
        else:
            GM_BGgenes.append(dist_df.loc[gene,case_compare])
            
    
    # print(len(GM_BGgenes), len(GM_genes))
    if stat_meas == 'MWU':
        statmwu,pvalmwu = mannwhitneyu(GM_genes, GM_BGgenes, alternative='greater')
    elif stat_meas == 'KST':
        statmwu,pvalmwu = kstest(GM_genes, GM_BGgenes, alternative='greater')

    fold_dif = np.median(GM_genes)/np.median(GM_BGgenes)
    print(fold_dif)

    fig, ax = plt.subplots(figsize=(11, 9))
    plt.hist(GM_BGgenes, bins = 100, label=f'(1) BG ({len(GM_BGgenes)})', alpha = 0.5, color='cornflowerblue')#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 
    plt.hist(GM_genes, bins = 100, label=f'(2) {PD_Gs} ({len(GM_genes)})', alpha = 0.5, color='red')#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 

    ax.set_yscale('log')
    title =  f'D{case_compare}'

    plt.title(title, fontsize=26) 
    ax.set_xlabel('Gene movement', fontsize=26) #dist_measure_l[dist_meas]
    ax.set_ylabel('frequency', fontsize=26)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=26)
        
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu <= 0.05:
        if stat_meas == 'MWU':
            fig.text(0.5, 0.52, f'MWU (2) > (1)\np-value = {pvalmwu_s}*\n', fontsize=28)
        elif stat_meas == 'KST':
            fig.text(0.5, 0.52, f'KST (2) > (1)\np-value = {pvalmwu_s}*\n', fontsize=28)
    else:
        if stat_meas == 'MWU':
            fig.text(0.5, 0.52, f'MWU (2) > (1)\np-value = {pvalmwu_s}\n', fontsize=28)
        elif stat_meas == 'KST':
            fig.text(0.5, 0.52, f'KST (2) > (1)\np-value = {pvalmwu_s}\n', fontsize=28)
    
    plt.savefig(f'{save_dir}/{case_compare}_{PD_Gs}vsOG_{stat_meas}.jpg', dpi = 350, format='jpg')                        
    plt.show()
    plt.close()
    
             

def plotPosDist(array_GM, title, save_dir, filename, S=1, v_line=None):  
    order = []
    array_GM = list(array_GM)
    array_GM.sort(reverse=True)
    array_GM = np.array(array_GM)
                    
    for i, gm in enumerate(array_GM):
        order.append(i) 
    
    
    fig, ax = plt.subplots(figsize=(8, 7))       
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    
    order_log = [i+1 for i in order]

    ax.set_xscale('log')

    ax.set_xlabel('Ranking', fontsize = 24)
    ax.set_ylabel('Total gene movement (TGM)', fontsize = 24)
    if v_line != None:
        plt.axvline(x = v_line, color = 'r',linewidth=6)
        fig.text(0.57,0.7,v_line,fontsize = 30)

    else:
        kneedle = KneeLocator(order, array_GM, S=S, curve="convex", direction="decreasing")
        elbowp = round(kneedle.elbow, 3)
        yval = array_GM[elbowp]
        
    colors = ['grey' if val < yval else default_colors[0] for val in array_GM]
    
    plt.scatter(order_log, array_GM, marker='o',#markerp[group[i]]
                          c=colors, 
                          alpha = 0.8,
                          s = 70)
        
    plt.axhline(y = yval, color = 'r',linewidth=3)

    fig.text(0.52,0.82,f'TGM threshold = {round(yval,3)}',fontsize = 18)
    fig.text(0.52,0.76,f'Nº predictions = {elbowp}',fontsize = 18)

    TP = elbowp/len(genes)*100
    sd2 = f'{save_dir}/{round(TP,3)}perc'
    if not os.path.exists(sd2):
        os.makedirs(sd2)  

    plt.savefig(f'{sd2}/{filename}',  dpi = 350, format='jpg', bbox_inches='tight')  
    plt.show()
    plt.close() 
    
    if v_line == None:
        return elbowp, yval

def euclidean_distance(df1, df2):
    distances = {}
    for column in df1.columns:
        distances[column] = [distance.euclidean(df1[column], df2[column])]
    return distances
   


# Normalized Euclidean distance - not called taht, use instead eculdean distance normalized by dividing with no
sensitivity = float(sys.argv[1])
Flag = str(sys.argv[2])
Flag2 = str(sys.argv[3])
wd = str(sys.argv[4])


with open(f'{wd}/input/PDgenes_DGN_ALL.pkl', 'rb') as handle:
    DGN_PD = pickle.load(handle) 


in_dir = 'step6_ComputeGMMs/output'
day_sets = [[8,18,25,32,37]]


if Flag:
    for days in day_sets:
        print(days)
        # days= [8, 18, 25, 32, 37]
        GM_stages = []    
        save_ext = '_'.join(f'D{str(x)}' for x in days)
    
        sd = f'{wd}/output/{save_ext}'
        if not os.path.exists(sd):
            os.makedirs(sd)  
       
        All_genes = []
        for file in os.listdir(f'{wd}/input/Genelists'):
            day = file.split('_')[2][:-4]
            if int(day) in days:
                Geneslist = pd.read_csv(f'{wd}/input/Genelists/{file}')
                Geneslist = Geneslist['genes'].tolist()
                All_genes.append(Geneslist)
        All_genes_inters = list(set.intersection(*map(set,All_genes)))
        All_genes_union = list(set.union(*map(set,All_genes)))
        pdgenes_dgn = list(set(DGN_PD)&set(All_genes_union))
        print(len(pdgenes_dgn))
        
        for day in days: 
            print(day)
    
            case_compare = f'D{day}'
            EucDist_PD = np.load(f'{in_dir}/ND_{day}_PWEucDist.npy')
            Genes_PD = pd.read_csv(f'{wd}/input/Genelists/Geneslist_ND_{day}.csv', header = 0)['genes'].values
            EucDist_PD_df = pd.DataFrame(EucDist_PD, index = Genes_PD, columns = Genes_PD)
            
            EucDist_C = np.load(f'{in_dir}/WT_{day}_PWEucDist.npy')
            Genes_C = pd.read_csv(f'{wd}/input/Genelists/Geneslist_WT_{day}.csv', header = 0)['genes'].values
            EucDist_C_df = pd.DataFrame(EucDist_C, index = Genes_C, columns = Genes_C)
            
            EucDist_PD_df = EucDist_PD_df[All_genes_inters]
            EucDist_PD_df = EucDist_PD_df.loc[All_genes_inters]
            EucDist_C_df = EucDist_C_df[All_genes_inters]
            EucDist_C_df = EucDist_C_df.loc[All_genes_inters]   
            
            dists = euclidean_distance(EucDist_PD_df, EucDist_C_df)
            GM = pd.DataFrame.from_dict(dists, orient='index',columns=[day])
            GM = GM.sort_values(by=day, ascending=False)
            
            GM_stages.append(GM)
          
            # PDGs2Background_GM(GM, case_compare, PDmap_noUK, sd, sim_measure = sim_meas, dist_meas = dist_meas, stat_meas = 'MWU', PD_Gs = 'PDmap_noUK')
            # PDGs2Background_GM(GM, case_compare, PDmap_PINK1, sd, sim_measure = sim_meas, dist_meas = dist_meas, stat_meas = 'MWU', PD_Gs = 'PDmap_PINK1')
           #  PDGs2Background_GM(GM, day, PDmap, sd, sim_measure = sim_meas, dist_meas = dist_meas, stat_meas = 'MWU', PD_Gs = 'PDmap')
            PDGs2Background_GM(GM, day, DGN_PD, sd, stat_meas = 'MWU', PD_Gs = 'DGN_PD')
    
        # GM across a set of stages              
        GM_stages_df = GM_stages[0].loc[All_genes_inters]
        for df in GM_stages[1:]:
            GM_stages_df = pd.merge(GM_stages_df, df.loc[All_genes_inters], left_index=True, right_index=True)
        print(GM_stages_df)
        GM_stages_df.to_csv(f'{sd}/GM_{save_ext}.csv')


          
     
if Flag2:  
    for days in day_sets:   
        # days = [18,25,32,37]
        print(days)
        save_ext = '_'.join(f'D{str(x)}' for x in days)
        sd = f'{wd}/output/{save_ext}'

        GM_stages_df = pd.read_csv(f'{sd}/GM_{save_ext}.csv',index_col=0,header=0)

        GM_stages_a = np.array(GM_stages_df.values)
        genes = GM_stages_df.index
        

    
        # Lenght of the GM vector in the stage-dim space
        GM_length = np.linalg.norm(GM_stages_a, axis=1)
        GM_length_df =  pd.DataFrame(GM_length, index = genes, columns=['GM_length']).sort_values(['GM_length'], ascending=False)                          
        genes = list(GM_length_df.index)


        
        elbow, yval = plotPosDist(GM_length, 'Thresholding genes according to TGM', sd, f'GM_length_{save_ext}_pos.jpg', sensitivity)
        TP = elbow/len(genes)*100
        print(TP)
        
 
        sd2 = f'{sd}/{round(TP,3)}perc'
        if not os.path.exists(sd2):
            os.makedirs(sd2)  
        TP_len = elbow
        print(TP_len)
        
        # topTP perc of genes based on GM length 
        TP_PDpreds = genes[:TP_len]          
        v_line = GM_length_df.loc[TP_PDpreds[TP_len-1],'GM_length']
        plot_hist_GML(GM_length, v_line, TP, TP_len, save_ext, sd2, save_ext)  
        # plotPosDist(GM_length, 'Gene movement length', sd2, f'GM_length_{save_ext}_pos.jpg', sens, TP_len)

        ribos_genes = []
        for gene in TP_PDpreds:
            if 'RPL' in gene or 'RPS' in gene:
                ribos_genes.append(gene)
        print(100*len(ribos_genes)/len(TP_PDpreds))
            
        # save dist of SS_PDpreds  
        GM_stages_PDpreds_df = GM_stages_df.loc[TP_PDpreds]
        GM_stages_PDpreds_df.to_csv(f'{sd2}/GM_{save_ext}_{round(TP,3)}_PDpreds.csv')
        GM_length_TP = GM_length_df.loc[TP_PDpreds].sort_values(['GM_length'], ascending=False)
        GM_length_TP.to_csv(f'{sd2}/GM_Length_{save_ext}_{round(TP,3)}_PDpreds.csv', index = True, header = True)
              
os.chdir(wd)                   
              
                
              
       
                
              