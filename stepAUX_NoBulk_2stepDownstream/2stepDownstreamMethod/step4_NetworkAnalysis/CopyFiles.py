# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 14:12:56 2023

@author: kmihajlo
"""

from shutil import copyfile
import os

work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')


copy_from = 'step2_StageSpecPreds/output'
copy_to = f'{work_dir}/input'  
file = 'CorePreds.pkl'
copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')



path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

copy_from = 'step2_NMTF_WeightedNets/output'
copy_to = f'{work_dir}/input/Genelists'                
files =  os.listdir(copy_from)    
for file in files:
    for file in files:
        if 'Geneslist_' in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')


# copy_from = 'step1_CreateNetworks/input/MolecularNetworks'
# copy_to = f'{work_dir}/input/MolecularNetworks' 
# if not os.path.exists(copy_to):
#     os.makedirs(copy_to) 
copy_to = f'{work_dir}/input'     
# copyfile('step11_ComparisonDEGsDAPs/output/D8_D18_D25_D32_D37/DAPs_D8_D18_D25_D32_D37_1FC.csv', f'{copy_to}/DAPs_D8_D18_D25_D32_D37_1FC.csv')
# copyfile('step11_ComparisonDEGsDAPs/output/D8_D18_D25_D32_D37/DEGs_D8_D18_D25_D32_D37_0.5FC.csv', f'{copy_to}/DEGs_D8_D18_D25_D32_D37_0.5FC.csv')
  
copy_from = 'step7_TGM_Predictions/output/D8_D18_D25_D32_D37/1.533perc'
file = 'GM_D8_D18_D25_D32_D37_1.533_PDpreds.csv'
copyfile(f'{copy_from}/{file}', f'{copy_to}/MO_{file}')

      
# copyfile(f'{copy_from}/Human_COEX_General_weight.edgelist', f'{copy_to}/Human_COEX_General_weight.edgelist')
# copyfile(f'{copy_from}/Human_GI_General.edgelist', f'{copy_to}/Human_GI_General.edgelist')
# copyfile(f'{copy_from}/Human_PPI_General.edgelist', f'{copy_to}/Human_PPI_General.edgelist')
# copyfile('step1_CreateNetworks/input/Metabolomics/enzyme_subprodgene_KEGG.csv', f'{copy_to}/enzyme_subprodgene_KEGG.csv')
# copyfile('step1_CreateNetworks/output/Entrez2Sym.pkl', f'{work_dir}/input/Entrez2Sym.pkl')
copyfile('step12_NetworkAnalysis/output/PPI_General_Biogrid_GeneSym.edgelist', f'{copy_to}/PPI_General_Biogrid_GeneSym.edgelist')



# copy_from = 'Data/ParsingGeneAnnotDBs/output'
# copyfile(f'{copy_from}/HSA_GO-BP.lst', f'{copy_to}/HSA_GO-BP.lst')
# copyfile(f'{copy_from}/HSA_Kegg_Pathways.lst', f'{copy_to}/HSA_Kegg_Pathways.lst')
# copyfile(f'{copy_from}/HSA_Reactome_Pathways.lst', f'{copy_to}/HSA_Reactome_Pathways.lst')
# copyfile(f'{copy_from}/HSA_KEGG_Pathways_meaning.lst', f'{copy_to}/HSA_KEGG_Pathways_meaning.lst')
# copyfile(f'{copy_from}/HSA_Reactome_Pathways_meaning.lst', f'{copy_to}/HSA_Reactome_Pathways_meaning.lst')
# copyfile(f'{copy_from}/Reactome_Group.csv', f'{copy_to}/Reactome_Group.csv')


# copyfile('Data/ParsingGeneAnnotDBs/go-basic.obo', f'{copy_to}/go-basic.obo')
# copyfile('step1_CreateNetworks/output/Entrez2Sym.pkl', f'{copy_to}/Entrez2Sym.pkl')
# copyfile('Data/PDmap/output/PDmap_BasePathways_noSBUK.lst', f'{work_dir}/input/PDmap_BasePathways_noSBUK.lst')             


# copyfile('step9_PubMedDEG_Validation/output/CommonlyExpressedGenes.txt', f'{copy_to}/CommonlyExpressedGenes.txt')
# copyfile('step10_KEGGEnrichmentAnalysis/output/All_genes_Union.txt', f'{copy_to}/All_genes_Union.txt')
# copyfile('CommonData/DEGs_Skupin.txt', f'{copy_to}/DEGs_Skupin.txt')
# copyfile('CommonData/LitValid_AllGenes.pkl', f'{copy_to}/LitValid_AllGenes.pkl')

os.chdir(work_dir) 
