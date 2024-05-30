# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 11:00:57 2021

@author: kmihajlo
"""
import pickle

PD_Genes = []

#parse t_cnv20200820
filename = 'input/t_cnv20200820.txt'
with open(filename, 'r', newline='\r\n') as f:
    next(f)
    for line in f:
        lspt = line.split('\t')
        PD_gene = lspt[0]
        PD_Genes.append(PD_gene)
        
#parse t_common_variant20200820
filename = 'input/t_common_variant20200820.txt'
with open(filename, 'r') as f:
    next(f)
    for line in f:
        lspt = line.split('\t')
        PD_gene = lspt[1]
        if ',' in PD_gene:
            PD_genes_many = PD_gene.split(',')
            for gene in PD_genes_many:
                PD_Genes.append(gene)
        if ';' in PD_gene:
            PD_genes_many = PD_gene.split(';')
            print(PD_genes_many)
            for gene in PD_genes_many:
                PD_Genes.append(gene)
        else:
            PD_Genes.append(PD_gene)
        
#parse t_dna_methylation20200820
filename = 'input/t_dna_methylation20200820.txt'
with open(filename, 'r', newline='\r\n') as f:
    next(f)
    for line in f:
        lspt = line.split('\t')
        PD_gene = lspt[3]
        if PD_gene != '':
            PD_Genes.append(PD_gene)
        
        
#parse t_gene_expression20200820
filename = 'input/t_gene_expression20200820.txt'
with open(filename, 'r', newline='\r\n') as f:
    next(f)
    for line in f:
        lspt = line.split('\t')
        PD_gene = lspt[1]
        PD_Genes.append(PD_gene)

#parse t_rare_gene20200820
filename = 'input/t_rare_gene20200820.txt'
with open(filename, 'r', newline='\r\n') as f:
    next(f)
    for line in f:
        lspt = line.split('\t')
        PD_gene = lspt[0]
        PD_Genes.append(PD_gene)
        
#parse t_rare_variant20200820
filename = 'input/t_rare_variant20200820.txt'
with open(filename, 'r', newline='\r\n') as f:
    next(f)
    for line in f:
        lspt = line.split('\t')
        #print(lspt)
        PD_gene = lspt[9]
        if ';' in PD_gene:
            PD_genes_many = PD_gene.split(';')
            PD_Genes.append(PD_genes_many[0])
        else:
            PD_Genes.append(PD_gene)
        
PD_Genes = list(set(PD_Genes))
with open('output/Gene4PD.pkl', 'wb') as handle:
	pickle.dump(PD_Genes, handle)  
    
    
    
    
    
    
    
    
    
    
    
    
    