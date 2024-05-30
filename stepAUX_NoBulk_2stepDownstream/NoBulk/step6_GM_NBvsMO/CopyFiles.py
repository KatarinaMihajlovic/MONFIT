# -*- coding: utf-8 -*-
"""
Created on Thu Oct  14 12:48:56 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os

work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')
      

copy_from = 'step1_NMTF/output'
copy_to = f'{work_dir}/input/Genelists'                
files =  os.listdir(copy_from)    
for file in files:
    for file in files:
        if 'Geneslist_' in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')

copy_to = f'{work_dir}/input'     

perc = '1.094'
copy_from = 'step4_TGM_Predictions/output/D8_D18_D37'
copyfile(f'{copy_from}/GM_D8_D18_D37.csv', f'{copy_to}/NB_GM_D8_D18_D37.csv')
copyfile(f'{copy_from}/{perc}perc/GM_D8_D18_D37_{perc}_PDpreds.csv', f'{copy_to}/NB_GM_D8_D18_D37_{perc}_PDpreds.csv')


path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

copy_from = 'step7_TGM_Predictions/output/D8_D18_D25_D32_D37'
copyfile(f'{copy_from}/GM_D8_D18_D25_D32_D37.csv', f'{copy_to}/MO_GM_D8_D18_D25_D32_D37.csv')
file = 'GM_D8_D18_D25_D32_D37_1.533_PDpreds.csv'
copyfile(f'{copy_from}/1.533perc/{file}', f'{copy_to}/MO_{file}')
file = 'GM_Length_D8_D18_D25_D32_D37_1.533_PDpreds.csv'
copyfile(f'{copy_from}/1.533perc/{file}', f'{copy_to}/MO_{file}')

# copyfile(f'{copy_from}/GM_D8_D18_D25_D32_D37.csv', f'{copy_to}/MO_GM_D8_D18_D25_D32_D37.csv')
# copyfile(f'{copy_from}/1.157perc/GM_D8_D18_D25_D32_D37_1.157_PDpreds.csv', f'{copy_to}/MO_GM_D8_D18_D25_D32_D37_1.157_PDpreds.csv')

           
copy_from = 'Data/ParsingGeneAnnotDBs/output'
copyfile(f'{copy_from}/HSA_GO-BP.lst', f'{copy_to}/HSA_GO-BP.lst')
copyfile(f'{copy_from}/HSA_Kegg_Pathways.lst', f'{copy_to}/HSA_Kegg_Pathways.lst')
copyfile(f'{copy_from}/HSA_Reactome_Pathways.lst', f'{copy_to}/HSA_Reactome_Pathways.lst')
copyfile(f'{copy_from}/HSA_KEGG_Pathways_meaning.lst', f'{copy_to}/HSA_KEGG_Pathways_meaning.lst')
copyfile(f'{copy_from}/HSA_Reactome_Pathways_meaning.lst', f'{copy_to}/HSA_Reactome_Pathways_meaning.lst')
copyfile(f'{copy_from}/Reactome_Group.csv', f'{copy_to}/Reactome_Group.csv')


copyfile('Data/ParsingGeneAnnotDBs/input/go-basic.obo', f'{copy_to}/go-basic.obo')
copyfile('step1_CreateNetworks/output/Entrez2Sym.pkl', f'{copy_to}/Entrez2Sym.pkl')
# copyfile('Data/PDmap/output/PDgenes_PDmap_noUK.pkl', f'{work_dir}/input/PDgenes_PDmap_noUK.pkl') 
copyfile('Data/PDmap/output/PDmap_BasePathways_noSBUK.lst', f'{work_dir}/input/PDmap_BasePathways_noSBUK.lst') 
# copyfile('Data/ParkinsonsUKAnnotationInitiative/output/PD_GOterms_human_UCL.lst', f'{work_dir}/input/PD_GO_UCL.lst') 
# copyfile('Data/PDmap/output/PDgenes_PDmap.pkl', f'{work_dir}/input/PDgenes_PDmap.pkl') 
copyfile('Data/DisGeNet/output/PDgenes_DGN_ALL.pkl', f'{work_dir}/input/PDgenes_DGN_ALL.pkl') 

# copyfile('Data/CorePredictions_prevStudy/CorePreds_GOtermsMeanining.txt', f'{work_dir}/input/CorePreds_GOtermsMeanining.txt')
# copyfile('Data/CorePredictions_prevStudy/CorePreds_KPmeaning.txt', f'{work_dir}/input/CorePreds_KPmeaning.txt')
# copyfile('Data/CorePredictions_prevStudy/CorePreds_RPtermsMeanining.txt', f'{work_dir}/input/CorePreds_RPtermsMeanining.txt')
# copyfile('PD_genesGOenrich/output/PDgenes_KPmeaning.txt', f'{copy_to}/PDgenes_KPmeaning.txt')
# copyfile('PD_genesGOenrich/input/PD_genes_DGN_DEGs.pkl', f'{copy_to}/PD_genes_DGN_DEGs.pkl')

os.chdir(work_dir)
