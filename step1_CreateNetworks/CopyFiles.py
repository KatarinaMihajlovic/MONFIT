from shutil import copyfile, copytree, rmtree
import os

work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 



if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')


copy_dir = 'Data/PD_Multiomics_ND_cellline/Transcriptomics'
copy_to = f'{work_dir}/input/Transcriptomics'

if not os.path.exists(copy_to):
    os.makedirs(copy_to)
    
if os.path.exists(copy_dir): 
    rmtree(copy_to) 
copytree(copy_dir, copy_to)

copy_from= 'Data/PD_Multiomics_ND_cellline/Proteomics/output'
copy_to = f'{work_dir}/input/Proteomics'
if not os.path.exists(copy_to):
    os.makedirs(copy_to)
copyfile(f'{copy_from}/Proteomicslog.csv', f'{copy_to}/Proteomicslog.csv')
    

copy_to = f'{work_dir}/input/Metabolomics'
if not os.path.exists(copy_to):
    os.makedirs(copy_to)
 

copy_from = 'Data/PD_Multiomics_ND_cellline/Metabolomics/output'    
copyfile(f'{copy_from}/Metabolomics_merged_KEGGid.csv', f'{copy_to}/Metabolomics_merged_KEGGid.csv')
copyfile(f'{copy_from}/KEGGID2metab.pkl', f'{copy_to}/KEGGID2metab.pkl')


copyfile('Data/Homo_sapiens.gene_info', f'{work_dir}/input/Homo_sapiens.gene_info')

copy_from = 'Data/MolecularNetworks'
copy_to = f'{work_dir}/input/MolecularNetworks'
if not os.path.exists(copy_to):
    os.makedirs(copy_to)
    
copyfile(f'{copy_from}/output/Human_PPI_General.edgelist', f'{copy_to}/Human_PPI_General.edgelist')
copyfile(f'{copy_from}/output/Human_COEX_General.edgelist', f'{copy_to}/Human_COEX_General.edgelist')
copyfile(f'{copy_from}/output/Human_GI_General.edgelist', f'{copy_to}/Human_GI_General.edgelist')
copyfile(f'{copy_from}/output/enzyme_subprodgene_KEGG.csv', f'{copy_to}/enzyme_subprodgene_KEGG.csv')

os.chdir(work_dir) 
