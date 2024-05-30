# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:41:58 2022

@author: kmihajlo
"""
import pandas as pd
import os, requests, re, pickle
from Bio import Entrez
import tqdm

def Genesym2Ensembl(filename = 'input/Homo_sapiens.gene_info'):
    Genesym2Ensembl = {}
    with open(filename) as f:
        f.readline()
        for line in f:
            lspt = line.split('\t')
            GeneNames = lspt[5].split('|')
            for Genename in GeneNames:
                if 'Ensembl' in Genename:
                    #print(Genename)
                    Ensembl = Genename.split(':')[1]
                    Genesym2Ensembl[lspt[2]] = Ensembl
    return Genesym2Ensembl

def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)

Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)
Genesym2Ensembl = Genesym2Ensembl()


def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='100000',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results   

def query_gene(gene, search_term):    
    counts = 0
    query  = f'({gene}[Text Word]) AND ({search_term}[Text Word])'
    query_file    = search(query)
    citation_list = query_file["IdList"]    
    counts = len(citation_list)   
    return counts





All_genes = []
for file in os.listdir('input/Genelists'):
    Geneslist = pd.read_csv(f'input/Genelists/{file}')
    Geneslist = Geneslist['genes'].tolist()
    All_genes.append(Geneslist)
# All_genes = set.intersection(*map(set,All_genes))

All_genes = set().union(*All_genes)
All_genes = set(All_genes)



LitValid_AllGenes = {} 
i = 0
for gene in tqdm.tqdm(All_genes):
    if gene not in LitValid_AllGenes:
        numPublications1 = query_gene(gene, 'Parkinson\'s disease') #instead of exctract_and
        try:
            numPublications2 = query_gene(Genesym2Ensembl[gene], 'Parkinson\'s disease') #instead of exctract_and
        except KeyError:
            numPublications2 = 0            
        numPublications = numPublications1 + numPublications2
        print(i, gene, numPublications)

        LitValid_AllGenes[gene] = [numPublications]
    i+=1

with open('output/LitValid_AllGenes.pkl', 'wb') as fp:   
    pickle.dump(LitValid_AllGenes, fp)
    
    