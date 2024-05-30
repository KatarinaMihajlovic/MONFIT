# -*- coding: utf-8 -*-
"""
.Created on Fri Oct 14 11:21:18 2022

@author: kimihajlo
"""

import csv, pickle, os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from sqlite3 import Error
import sqlite3
import pandas as pd


def write_txt(list_l, file_path, name):
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    with open(f'{file_path}/{name}.txt', 'w') as f: 
        for x in list_l:
            f.write(f'{x}\n')
            
def createConnection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn

def showAllTables(conn):
    """
    Query all tabels
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute('SELECT name from sqlite_master where type= "table"')

    rows = cur.fetchall()

    return rows

def selectAllFromTable(conn, t):
    """
    Query all rows in the tasks table
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute("SELECT * FROM {}".format(t))

    rows = cur.fetchall()

    return rows

def selectHeaderTable(conn, t):

    cur = conn.cursor()
    cur.execute("PRAGMA table_info({})".format(t))

    rows = cur.fetchall()
    return rows


def parse_DisGeNet_summary(file_name, threshold):
    PD_DisGeNetgenes_summary = open(file_name, encoding = 'utf-8')
    read_tsv = csv.reader(PD_DisGeNetgenes_summary, delimiter="\t")
    next(read_tsv, None)
    
    PD_DisGeNetgenes = []
    PD_DisGeNet_onlygenes_symb = []
    PD_DisGeNet_onlygenes_Entrez = []
    PD_gene_scores = []

    
    for line in read_tsv:
        PD_gene = line[2]
        PD_gene_Ent = line[3]
        PD_gene_score = line[11]
        PD_gene_mentions = line[14]
        
        if float(PD_gene_score) > threshold:
            PD_gene_scores.append(PD_gene_score)
            PD_DisGeNet_onlygenes_symb.append(PD_gene)
            PD_DisGeNet_onlygenes_Entrez.append(PD_gene_Ent)
            PD_DisGeNetgenes.append([PD_gene, PD_gene_score, PD_gene_mentions])

    return PD_DisGeNetgenes

def parse_DisGeNet_evidence(file_name, threshold):
    cats_dict = {'CausalMutation':[], 'Therapeutic':[], 
                 'Biomarker':[], 'GeneticVariation':[]}
    
    PD_DisGeNetgenes_evidences = open(file_name, encoding = 'utf-8')
    read_tsv = csv.reader(PD_DisGeNetgenes_evidences, delimiter="\t")
    next(read_tsv, None)
    
    PD_gene_cats = []
    PDgenes_DB = []
    PD_genes = []
    for line in read_tsv:
        PD_gene = line[2]
        PD_gene_score = line[4]
        PD_gene_cat = line[5]
        PD_gene_cats.append(PD_gene_cat)
        PDgene_DB = line[7]

        if float(PD_gene_score) > threshold and PDgene_DB != 'BEFREE' and PDgene_DB != 'LHGDN':
            cats_dict[PD_gene_cat].append(PD_gene)
            PDgenes_DB.append(PDgene_DB)
            PD_genes.append(PD_gene)
            
    for cat in cats_dict:
        cats_dict[cat] = list(set(cats_dict[cat]))
    # print(set(PD_gene_cats))
    PDgenes_DB_dict = dict(zip(PD_genes, PDgenes_DB))

    return cats_dict, PDgenes_DB_dict



#### Query DGN database
Flag = True
if Flag == True:
    conn = createConnection('input/disgenet_2020.db')
    
    t1 = 'geneDiseaseNetwork'
    t2 = 'geneAttributes'
    t3 = 'diseaseAttributes'
    t4 = 'disease2class'
    t5 = 'diseaseClass'
    join12 = 'geneNID'
    join13 = 'diseaseNID'
    join34 = 'diseaseNID'
    join45 = 'diseaseClassNID'
    
    q = "SELECT * FROM (SELECT * FROM {t1} LEFT JOIN {t2} ON {t1}.{join12} = {t2}.{join12}) AS T12  LEFT JOIN (SELECT * FROM {t3} LEFT JOIN (SELECT *FROM {t4} LEFT JOIN {t5}  ON {t4}.{join45} = {t5}.{join45} ) AS T45 ON {t3}.{join34} = T45.{join34} ) AS T345 ON T12.{join13} = T345.{join13}".format(t1=t1, t2=t2, t3=t3, t4=t4, t5=t5, join12=join12, join13=join13, join34=join34, join45=join45)
    
    
    cur = conn.cursor()
    cur.execute(q)
    
    rows = cur.fetchall()
    
    headerAll = []
    for t in [t1,t2,t3,t4,t5]:
        headerAll.extend(list(pd.DataFrame(selectHeaderTable(conn, t))[1].unique()))
    
    
    dfAll = pd.DataFrame(rows, columns = headerAll)
    
    dfAll.columns = headerAll
    dfAll = dfAll.loc[:,~dfAll.columns.duplicated()]
    dfAll.to_csv('output/SQLDB_gene_disease_associations.csv', index=False)
            
 
# Get all genes related to PD   
DGN_DB = pd.read_csv('output/SQLDB_gene_disease_associations.csv', header=0, low_memory=False)

#save all DGN PD genes
genes = list(set(DGN_DB.loc[DGN_DB['diseaseId'] == 'C0030567', 'geneName'].tolist()))
with open('output/PDgenes_DGN_ALL.pkl', 'wb') as handle:
    pickle.dump(genes, handle)  







        
        
        