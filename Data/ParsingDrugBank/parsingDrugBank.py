import xml.etree.ElementTree as ET
import pandas as pd


def parsingXML2Targets(filePath, savePath = ''):
    """
    Pasing DrugBank XML to CSV with Targets of each Drugs
    
    Parameter
    ---------
    filePath: str
        Path with the xml file
        
    Return
    ------
        DataFrame with the drugs and the corresponding targets
    """
    
    tree = ET.parse(filePath)
    root = tree.getroot()

    res = []
    for drugs in root:
        t = drugs.attrib['type']
        for atr in drugs:
           
            if ('drugbank-id' in atr.tag)&('primary' in atr.attrib):
                dbID = atr.text
                
            if ('name' in atr.tag):
                dbName = atr.text
    
            if ('groups' in atr.tag):
                for groups in atr:
                    if 'group' in groups.tag:
                        group = groups.text   
    
            if ('targets' in atr.tag):
                for targets in atr:
                    for tInfo in targets:
    
                        if 'polypeptide' in tInfo.tag: 
                            geneTarget = ''
                            for tInfoGene in tInfo:
                                if ('gene-name' in tInfoGene.tag)&(tInfoGene.text is not None): 
                                    geneTarget = tInfoGene.text
                                if ('organism' in tInfoGene.tag)&(geneTarget != ''):
                                    if("Humans" in tInfoGene.text):
                                        res.append([t, dbID, dbName, group, geneTarget])
    
    dfRes = pd.DataFrame(res, columns=['type', 'DrugBank ID', 'DrugBank Name', 'Group', 'Gene Target'])
    if savePath != '': dfRes.to_csv(savePath, index=False) 
    return dfRes



DTI = parsingXML2Targets('input/full database.xml', savePath = 'output/DTI.csv')
DTI_approved_invest = DTI[(DTI['Group'] == 'approved') | (DTI['Group'] == 'investigational') ]
DTI_approved_invest.to_csv('output/DTI_invest_approv.csv', index=False) 


