import os
import pandas as pd

#Parse NCBI_gene2go_GO-BP.txt -> gene - GO_term
GO_BP_input = open('input/gene2go','r')
GO_BP_input.readline()
GO_BP_output = open('output/HSA_GO-BP.lst','w')

for line in GO_BP_input.readlines():
    lspt = line.strip().split('\t')
    if lspt[0] == '9606' and lspt[7] == 'Process':
        GO_BP_output.write(lspt[1] + '\t' + lspt[2] + '\n')

GO_BP_input.close()
GO_BP_output.close()


#Parse NCBI2Reactome_All_Levels.txt -> gene - Reactome_pathwayID
REACTOME_input = open('input/NCBI2Reactome_All_Levels.txt','r')
# REACTOME_input.readline()
REACTOME_output = open('output/HSA_Reactome_Pathways.lst','w')
REACTOME_outputM = open('output/HSA_Reactome_Pathways_meaning.lst','w')

for line in REACTOME_input.readlines():
    lspt = line.strip().split('\t')
    # print(lspt)
    if 'HSA' in lspt[1]:
        REACTOME_output.write(lspt[0] + '\t' + lspt[1] + '\n')
        REACTOME_outputM.write(f'{lspt[1]}\t{lspt[3]}\n')
        
REACTOME_input.close()
REACTOME_output.close()


# Reactome hierarchy 
dfPathways = pd.read_csv('input/ReactomePathwaysRelation.txt', sep='\t', header=None)
dfPathwaysInfo = pd.read_csv('input/ReactomePathways.txt', sep='\t', names= ['Pathways', 'Name', 'Species'])

dfPathwaysHumanPathways = dfPathways[dfPathways[0].str.contains(r'\bR-HSA')]
parent, child = dfPathwaysHumanPathways[0].values, dfPathwaysHumanPathways[1].values
rootsPathways = set([p for p in parent if p not in child])

toSave = []
for pathway in rootsPathways:

    toSave.append([pathway, pathway, pathway])
    groups = dfPathways[dfPathways[0].isin([pathway])][1].values

    orderedListOfPathways = [pathway]

    for subgroup in groups:
        toSave.append([pathway, subgroup, subgroup])
        tmpPathways = [subgroup]
        newPathways = [subgroup]

        while (len(newPathways) != 0):
            newPathways = dfPathways[dfPathways[0].isin(newPathways)][1].values
            for p in newPathways:
                if p not in tmpPathways:
                    tmpPathways.append(p)
                    toSave.append([pathway, subgroup, p])
        orderedListOfPathways.extend(tmpPathways)

dfToSave = pd.DataFrame(toSave, columns=['Group', 'Subgroup', 'Pathway'])
dfTmp = pd.merge(dfToSave, dfPathwaysInfo, right_on='Pathways', left_on='Group')[['Name', 'Subgroup', 'Pathway']]
dfTmp.columns = ['Group', 'Subgroup', 'Pathway']
dfTmp = pd.merge(dfTmp, dfPathwaysInfo, right_on='Pathways', left_on='Subgroup')[['Group', 'Name', 'Pathway']]
dfTmp.columns = ['Group', 'Subgroup', 'Pathway']

dfTmp.to_csv('output/Reactome_Group.csv', index=False) 



#Parse KEGG_genes_pathways.txt -> gene - KEGG_pathwayID
KEGG_input = open('input/KEGG_path2gene.txt','r')
KEGG_output = open('output/HSA_KEGG_Pathways.lst','w')

for line in KEGG_input.readlines():
	lspt = line.strip().split('\t')
	gene = lspt[0].split(':')[1]
	pathway = lspt[1].split(':')[1]
	KEGG_output.write(pathway + '\t' + gene + '\n')

KEGG_input.close()
KEGG_output.close()


KEGG_input = open('input/hsa','r')
KEGG_output = open('output/HSA_KEGG_Pathways_meaning.lst','w')
for line in KEGG_input.readlines():
    lspt = line.strip().split('\t')
    pathway = lspt[0]
    name = lspt[1].replace(' - Homo sapiens (human)', '')
    KEGG_output.write(pathway + '\t' + name + '\n')

KEGG_input.close()
KEGG_output.close()