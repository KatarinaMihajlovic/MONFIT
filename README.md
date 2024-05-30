# MONFIT
Multi-Omics Non-negative matrix tri-Factorization Integration of Time-series
MONFIT uses NMTF to integrate multi-omics time-point-specific data of scRNA-seq, bulk proteomics, and bulk metabolomics measurements of disease and control samples from a time-series experiment with prior knowledge from molecular interaction networks in terms of PPI, metabolic-interaction (MI), genetic interaction (GI), and COEX, producing gene embeddings. Then, MONFIT holistically mines these embeddings across all time points to identify new disease-associated genes, by applying a new downstream method based on the concept introduced in Mihajlovic et al. [Link](https://doi.org/10.1038/s41598-024-61844-3), which states that those genes whose position changes the most between disease and control samples across all time points are disease-related. 

## What can you find in this repository?
This repository contains all data, scripts and results related to our recent work. 
In particular, you will find the Jupyter notebook [Tutorial.ipynb](https://github.com/KatarinaMihajlovic/MONFIT/blob/main/Tutorial.ipynb) contains the complete MONFIT pipeline and computational methods for analyzing the gene predictions.

### Additional Information
In "Data/MolecularNetworks" directory:

For constructing the PPI network, download BIOGRID-ALL-4.4.218.tab3.zip from BioGrid and unpack it in the "input" directory. Then, run ConstructPPI.py.

For constructing the GI network, download BIOGRID-ALL-4.4.218.tab3.zip from BioGrid and unpack in the "input" directory. Then, run ConstructGI.py.

For constructing the COEX network, download Hsa-u.v22-05.G16651-S245698.combat_pca.subagging.z.d.zip from Coexpresdb and unpack in "input/Hsa-u.v22-05.G16651-S245698.combat_pca.subagging.z.d" directory. Then, run ConstructCOEX.py.

For constructing the MI network in step1_CreateNetworks, download information on enzymes from KEGG using ParseKEGG_enzymeData.sh. We provide the list of human enzymes we downloaded from KEGG on 06.03.2023. Then, run MI_Enzymes2Metabolites2Genes.py to obtain the file connecting enzymes, metabolites and genes.

Note: To run the .sh script on Windows, it is necessary to use an environment that can interpret and execute the script (e.g., Git Bash). If using Git Bash, make sure that "wget" utility is installed on the system.

In the "output" directory, we provide PPI, GI, and COEX networks used in our study and KEGG data used for constructing the data-specific MI networks.

In "Data/parsingDrugBank" directory:
Download drugbank_all_full_database.xml file from go.drugbank.com; version 	5.1.10 released on 2023-01-04.

In "Data/ParsingGeneAnnotDBs" directory:
gene2go.gz downloaded from https://ftp.ncbi.nlm.nih.gov/gene/DATA/
NCBI2Reactome_All_Levels.txt, ReactomePathways.txt and ReactomePathwaysRelation.txt downloaded from https://reactome.org/download/current/
KEGG_PathwaysList.txt and KEGG_Pathway2Gene.txt downloaded using ParseKEGG.sh
go-basic.obo downloaded from http://geneontology.org/docs/download-ontology/
All files were downloaded on 21.03.2023.

In "Data/PD_Multiomics_ND_cellline/Transcriptomics" directory:
WT denotes Wild Type, i.e., control samples.
ND denotes samples of a PD cell line harbouring a PINK1 mutation.

In "step2_NMTF_WeightedNets" directory:
To ensure reproducibility, we provide matrix factors that the step 1 of MONFIT produces ("output"). Running NMTF on different systems may result in slight numerical differences.

In "step3_ClusterG1NMTF_WeightedNets" directory:
To ensure reproducibility, we provide k-means clusters, because the k-means clustering algorithm is non-deterministic.

In "stepAUX_NoBulk_2stepDownstream" directory:
We provide scripts for comparing the results obtained with MONFIT on scRNA-seq, bulk proteomics, and bulk metabolomics measurements of disease and control samples from a Parkinson's disease time-series experiment with prior knowledge from molecular interaction networks against:
- the predictions obtained by the downstream analysis approach that inspired MONFIT's pipeline, which is presented in Mihajlovic et al. [Link](https://doi.org/10.1038/s41598-024-61844-3).
- the predictions obtained by not including bulk metabolomics and proteomics data during MONFIT integration ("NoBulk" directory). Here, we also provide matrix factors that the step 1 of MONFIT produces.

In "stepAUX_RobustnessAnalysis_k1k2" and "stepAUX_RobustnessAnalysis_Weights" directories:
We provide scripts for assessing the robustness of MONFIT to the i) number of dimensions and ii) different combinations of weighting factors of the input matrices used for producing the lower-dimensional matrices with the step 1 of MONFIT.


### How to run the notebook
pip install -r requirements.txt

Execute the jupyter notebook [Tutorial.ipynb](https://github.com/KatarinaMihajlovic/MONFIT/blob/main/Tutorial.ipynb) 

