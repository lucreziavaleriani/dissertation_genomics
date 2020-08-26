#!usr/bin/python3
import sys
import shutil
import tempfile
import urllib.request
import ssl
import csv
import pandas as pd
import requests
from bs4 import BeautifulSoup as bs
ssl._create_default_https_context = ssl._create_unverified_context

#['#AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinicalSignificance', 'ClinSigSimple', 'LastEvaluated', 'RS#(dbSNP)', 'nsv/esv(dbVar)', 
# 'RCVaccession', 'PhenotypeIDS', 'PhenotypeList', 'Origin', 'OriginSimple', 'Assembly', 'ChromosomeAccession', 'Chromosome', 'Start', 'Stop', 'ReferenceAllele', 
# 'AlternateAllele', 'Cytogenetic', 'ReviewStatus', 'NumberSubmitters', 'Guidelines', 'TestedInGTR', 'OtherIDs', 'SubmitterCategories', 'VariationID'] #31 campi

def clean_clinvar(clinvar):
    cv = open(clinvar)
    orpha=[] #index della line che ha corrispondenza con orphanet
    omim=[]
    a=[]
    for line in cv:
        a.append(line)
        om=line.find('OMIM')
        orph=line.find('Orphanet')
        if orph != -1:
            orpha.append(line)
        if om != -1:
            omim.append(line)
    l = list(set(omim+orpha))
    return l

def filter_cv(cv_cl):
    colnames = ['#AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinicalSignificance', 'ClinSigSimple', 'LastEvaluated', 'SNP', 'nsv/esv(dbVar)', 'RCVaccession', 'PhenotypeIDS', 'PhenotypeList']
    data = pd.read_csv(cv_cl, delimiter=',', names=colnames)
    d={}
    c=0
    name = data.Name.tolist()
    ids = data.PhenotypeIDS.tolist()
    path=data.ClinicalSignificance.tolist()
    tp=data.Type.tolist()
    snp=data.SNP.tolist()
    for i in ids:
        c+=1
        n=name[c-1]
        if path[c-1]=='Pathogenic' and tp[c-1].strip('|')=='single|nucleotide|variant':  
            d[n]=d.get(n,{})
            d[n]['omim']=d[n].get('omim',[])
            d[n]['orphanet']=d[n].get('orphanet',[])
            d[n]['SNP']=d[n].get('SNP',[])
            d[n]['SNP'].append(snp[c-1])
            i=i.split(',')    
            for j in i:
                o=j.find('OMIM') 
                if o != -1:
                    om=j[o+5: o+11]
                    d[n]['omim'].append(om)
                r=j.find('Orphanet')
                if r != -1:
                    orph=j[r+14:]
                    f=orph.find('|')
                    if f!=-1:
                        orph=orph[:f]
                    d[n]['orphanet'].append(orph)
    return d



if __name__ == '__main__':
    # cv_cl=sys.argv[1]
    # cv_filter= filter_cv(cv_cl)
    clinvar = sys.argv[1]
    cv_d=clean_clinvar(clinvar)
    # with open('clinvar_clean.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in cv_d:
    #         l=i.replace(' ','|')
    #         l=l.replace(';','|')
    #         m=l.split()
    #         all_writer.writerow([m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11],m[12],m[13]])
    # with open('cv_snv_patho.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     c=','
    #     for i in cv_filter:
    #         all_writer.writerow([i,c.join(set(cv_filter[i]['omim'])),c.join(set(cv_filter[i]['orphanet'])),c.join(set(cv_filter[i]['SNP']))])