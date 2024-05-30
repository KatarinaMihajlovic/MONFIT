# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:44:51 2023

@author: kmihajlo
"""

import pyreadr, pickle
import pandas as pd
import numpy as np 

def avg_replicates(df, objects, cc):
    abundances = []
    for obj in objects:
        abundance_reps = np.array(df[obj])
        abundance_reps = np.append(abundance_reps,[0,0])
        zeros = len(abundance_reps) - np.count_nonzero(abundance_reps)
        
        if zeros == 0:
            abundance = np.mean(abundance_reps)
        elif zeros <= len(abundance_reps)/2:
            abundance_reps_no_0 = abundance_reps[abundance_reps!=0]
            abundance = np.mean(abundance_reps_no_0)
        else:
            abundance = 0
        abundances.append(abundance)
    
    df_avg = pd.DataFrame(abundances,columns=[f'{cc}'],index=objects)
    return df_avg

              
#metabolomics data        
result = pyreadr.read_r('input/M_GC_spread_WT_vs_ND_8to37.RData')      
Metabolomics_GC = result[None]        
result = pyreadr.read_r('input/M_LC_spread_WT_vs_ND_8to37.RData')      
Metabolomics_LC = result[None]        

# remove Gluconic acid 6-phosphate  form Metabolomics_GC, as it is the same as '6-Phosphogluconic acid', and L-Trypthophan is the same as Tryptophan
Metabolomics_GC = Metabolomics_GC.drop(['Gluconic acid 6-phosphate'], axis=1)
Metabolomics_LC = Metabolomics_LC.drop(['L-Tryptophan'], axis=1)

# print(Metabolomics_GC['Gluconic acid 6-phosphate'])
# rename metabolites, different technologies use different naming, so rename to have shared naming from GC
rename_map = {'a-Ketoglutarate':'2-Oxoglutaric acid', 'AMP':'Adenosine monophosphate', 'Aspartate':'Aspartic acid', 
              'Citrate':'Citric acid', 'Fumarate':'Fumaric acid', 'GABA':'Gamma-Aminobutyric acid', 'Glucose-6-phosphate':'Glucose 6-phosphate', 
              'Glutamate':'Glutamic acid', 'Lactate':'Lactic acid', 'N-acetylaspartate':'N-Acetyl-L-aspartic acid', 'UMP':'Uridine monophosphate'}

Metabolomics_GC = Metabolomics_GC.rename(columns=rename_map)
Metabolomics_LC = Metabolomics_LC.rename(columns=rename_map)

metab_GC = Metabolomics_GC.columns.to_list()[3:]        
metab_LC = Metabolomics_LC.columns.to_list()[3:]        
   
intersection = list(set(metab_GC)&set(metab_LC))     
# print(list(set(metab_GC)&set(metab_LC)))        

# scale each metabolomics dataset with the mean of non-zero values
metab_GC_vals = Metabolomics_GC[metab_GC].values        
metab_GC_no0_mean = metab_GC_vals[np.where(metab_GC_vals!=0)].mean() 
metab_GC_vals = metab_GC_vals/metab_GC_no0_mean        
Metabolomics_GC_scaled = Metabolomics_GC.copy()
Metabolomics_GC_scaled[metab_GC] = metab_GC_vals 

metab_LC_vals = Metabolomics_LC[metab_LC].values        
metab_LC_no0_mean = metab_LC_vals[np.where(metab_LC_vals!=0)].mean() 
metab_LC_vals = metab_LC_vals/metab_LC_no0_mean        
Metabolomics_LC_scaled = Metabolomics_LC.copy()
Metabolomics_LC_scaled[metab_LC] = metab_LC_vals 

days = set(Metabolomics_LC['Day'])
days = [int(x) for x in days]
types_cc = set(Metabolomics_LC['Condition'])

# average across replicates for each condition 
ccs = []
metabGC_scaled_avg = []
metabLC_scaled_avg = []

for type_cc in types_cc:
    for day in days:
        cc = f'{type_cc}_{day}'
        ccs.append(cc)
        cc_MGC = Metabolomics_GC_scaled.loc[(Metabolomics_GC_scaled['Condition'] == type_cc) & (Metabolomics_GC_scaled['Day'] == day)]
        cc_MLC = Metabolomics_LC_scaled.loc[(Metabolomics_LC_scaled['Condition'] == type_cc) & (Metabolomics_LC_scaled['Day'] == day)]
        
        cc_MGC_avg = avg_replicates(cc_MGC, metab_GC, cc)
        cc_MLC_avg = avg_replicates(cc_MLC, metab_LC, cc)
        
        metabGC_scaled_avg.append(cc_MGC_avg[cc].values)
        metabLC_scaled_avg.append(cc_MLC_avg[cc].values)
       
metabGC_scaled_avg_df = pd.DataFrame(metabGC_scaled_avg, index = ccs, columns = metab_GC)
metabLC_scaled_avg_df = pd.DataFrame(metabLC_scaled_avg, index = ccs, columns = metab_LC)

# drop duplicate metabolites  
metabGC_scaled_avg_df_unique = metabGC_scaled_avg_df.drop(intersection, axis=1)
metabLC_scaled_avg_df_unique = metabLC_scaled_avg_df.drop(intersection, axis=1)
 

# merge 2 metabolomics datasets
metab_final = pd.concat([metabGC_scaled_avg_df_unique, metabLC_scaled_avg_df_unique], axis=1)

# choose the values  for duplicate metabolites from 1 measuremnt type, based on the highest avg value across all cell conditions
metab_inter_vals = []
for metab in intersection:
    metabGC = list(metabGC_scaled_avg_df[metab].values)
    metabLC = list(metabLC_scaled_avg_df[metab].values)
    
    if np.mean(metabGC) > np.mean(metabLC):
        metab_inter_vals.append(metabGC)
    else:
        metab_inter_vals.append(metabLC)


metab_final[intersection] = np.array(metab_inter_vals).T
metab_final.to_csv('output/Metabolomics_merged.csv',index=True)  


mappingMets = pd.read_excel('input/MappingMetaboliteNames.xlsx')
mapDict_KEGG = dict(zip(mappingMets.Query,mappingMets.KEGG))

metab_KEGGmapped = metab_final.rename(columns = mapDict_KEGG)
metab_KEGGmapped = metab_KEGGmapped.loc[:, metab_KEGGmapped.columns.notna()]
metab_KEGGmapped.to_csv('output/Metabolomics_merged_KEGGid.csv',index=True)  


used_metabs = [x for x in metab_final.columns if not isinstance(mapDict_KEGG[x], float)]
finalKEGG_map = dict(zip(metab_KEGGmapped.columns,used_metabs))

with open('output/KEGGID2metab.pkl', 'wb') as f:
    pickle.dump(finalKEGG_map, f)


      
     
