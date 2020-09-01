#!usr/bin/python3
import sys
import shutil
import tempfile
import urllib.request
import ssl
import csv
import pandas as pd
import requests
from bs4 import  BeautifulSoup as bs
ssl._create_default_https_context = ssl._create_unverified_context


aa={'Lys':'K','Arg':'R','His':'H','Asp':'D','Glu':'E','Ala':'A','Val':'V','Ile':'I','Leu':'L',
'Met':'M','Phe':'F','Tyr':'Y','Trp':'W','Cys':'C','Sec':'U','Gly':'G','Pro':'P','Ser':'S','Thr':'T',
'Asn':'N','Gln':'Q'}

def omim_var(humvar):
    var = open(humvar)
    d_var = {}
    d_dis_acr = {}
    om=[]
    lvar=[]
    gene_up={}
    for line in var:
        if line[21:24] == 'VAR':
            gene = line[0:9].strip()
            AC = line[10:20].strip()
            gene_up[gene]=AC

        if line[21:24] == 'VAR' and line[77] != '-':
            gene = line[0:9].strip()
            AC = line[10:20].strip()
            VAR = line[21:32].strip()
            aa=line[33:47].strip()
            dbSNP = line[62:76].strip()
            disease = line[77:-1].strip()
            if disease[-11:-8] == 'MIM':
                b = disease.find(("("))
                dis = disease[:b-1]
                acr = disease[b+1:-14]
                d_dis_acr[dis] = acr
                mim = disease[-7:-1]
                om.append(mim)
                if [gene,aa] not in lvar:
                    lvar.append([gene,aa])
                d_var[acr] = d_var.get(acr, {})
                d_var[acr]['omim'] = d_var[acr].get('omim', mim)
                d_var[acr]['gene'] = d_var[acr].get('gene', gene)
                d_var[acr]['var'] = d_var[acr].get('var', [])
                d_var[acr]['var'].append(VAR)
                d_var[acr]['snp'] = d_var[acr].get('snp', [])
                d_var[acr]['snp'].append(dbSNP)
                d_var[acr]['change'] = d_var[acr].get('change', [])
                d_var[acr]['change'].append(aa)
    print(len(gene_up))
    return d_var,lvar,gene_up


def orphanet_omim(xml1):
    with open(xml1, "r", encoding='iso-8859-1') as file:
        all_orpha = {}
        all_omim = []
        file.readline()
        xml_soup = bs(file.read(), 'xml')
        result = xml_soup.find_all('Disorder')
        val = []
        for j in result:
            for i in j.find('OrphaNumber'):
                orpha = i.string
            result2 = j.find_all('ExternalReference')
            for i in result2:
                for a in i.find('Source'):
                    if a.string == 'OMIM':
                        for c in i.find('DisorderMappingValidationStatus'):
                            for d in i.find('Name'):
                                val = d.string
                                p = val.find('(')
                                val = val[:p-1]
                        for b in i.find('Reference'):
                            om = b.string
                            all_omim.append(om)
                            all_orpha[orpha] = all_orpha.get(
                                orpha, [[om], val])
                            if om not in all_orpha[orpha][0]:
                                all_orpha[orpha][0].append(om)
    print(len(all_orpha))
    print(len(set(all_omim)))
    return all_orpha


def orphanet_group(csv2):
    colnames = ['orphan', 'desc', 'xml']
    data = pd.read_csv(csv2, delimiter=';', names=colnames)
    xml_l = data.xml.tolist()
    orph = data.orphan.tolist()
    orph_list = orph[1:]
    xml_list = xml_l[1:]
    return orph_list, xml_list


def group_orphanet(orph_list, xml_list):
    group_orpha = {}
    al_orpha = []
    for i in range(len(orph_list)):
        x = xml_list[i]
        with urllib.request.urlopen(x) as response:
            with tempfile.NamedTemporaryFile() as tmp_file:
                shutil.copyfileobj(response, tmp_file)
                with open(tmp_file.name, "r", encoding='iso-8859-1') as file:
                    file.readline()
                    xml_soup = bs(file.read(), 'xml')
                    cla = xml_soup.find_all('ClassificationNode')
                    for j in cla:
                        if j.find('OrphaCode') != None :
                            for o in j.find('OrphaCode'):
                                orpha = o.string
                                al_orpha.append(orpha)
                                group_orpha[orph_list[i]] = group_orpha.get(
                                    orph_list[i], [])
                                group_orpha[orph_list[i]].append(orpha)
    # for i in group_orpha:
    #     print(i, len(set(group_orpha[i])))
    return group_orpha


def orphanet_var(mim_var, orphanet_mim):
    orph_var = {}
    #all_omim=[]
    all_var=[]
    dis=list(mim_var.keys())
    for i in orphanet_mim: #orpha
        for j in orphanet_mim[i][0]: #omim 
            for l in dis:
                if j == mim_var[l]['omim']:
                    orph_var[i] = orph_var.get(i, {})
                    
                    all_var.append([mim_var[l]['change'],mim_var[l]['gene'],i])

                    orph_var[i]['var']=orph_var[i].get('var',[])
                    if mim_var[l]['var'] not in orph_var[i]['var']:
                         orph_var[i]['var']+= mim_var[l]['var']
                   
                    orph_var[i]['snp']=orph_var[i].get('snp',[])
                    if mim_var[l]['snp'] not in orph_var[i]['snp']:
                         orph_var[i]['snp']+= mim_var[l]['snp']

                    orph_var[i]['omim']=orph_var[i].get('omim',[])
                    if j not in orph_var[i]['omim']:
                        orph_var[i]['omim'].append(j)

                    orph_var[i]['gene']=orph_var[i].get('gene',[])
                    if mim_var[l]['gene'] not in orph_var[i]['gene']:
                         orph_var[i]['gene'].append(mim_var[l]['gene']) 

                    orph_var[i]['change']=orph_var[i].get('change',[])
                    if mim_var[l]['change'] not in orph_var[i]['change']:
                         orph_var[i]['change']+= mim_var[l]['change']
                    
                    orph_var[i]['val']=orph_var[i].get('val',orphanet_mim[i][1])
    print(all_var)
    return orph_var,all_var

def group_var(orpha_var, group_orpha):
    g_var={}
    all_exact_var=[]
    for i in group_orpha:
        for j in orpha_var:
            if j in group_orpha[i]:
                g_var[i]= g_var.get(i,{})
                g_var[i]['var']=g_var[i].get('var',[])
                g_var[i]['var']+= orpha_var[j]['var']
                
                g_var[i]['snp']=g_var[i].get('snp',[])
                g_var[i]['snp']+= orpha_var[j]['snp']
                
                g_var[i]['orpha']=g_var[i].get('orpha',[])
                g_var[i]['orpha'].append(j)
                
                g_var[i]['gene']=g_var[i].get('gene',[])
                g_var[i]['gene']+= orpha_var[j]['gene']
                
                g_var[i]['change']=g_var[i].get('change',[])
                g_var[i]['change']+= orpha_var[j]['change']

                g_var[i]['val']=g_var[i].get('val',[])
                g_var[i]['val'].append(orpha_var[j]['val'])
                if orpha_var[j]['val'] == 'E':
                    all_exact_var+=orpha_var[j]['var']
    # for i in g_var:
    #     print(i, len(set(g_var[i]['var'])))
    return g_var

def orphanet_swiss(xml2):
    uniprot=[]
    with open(xml2, "r", encoding='iso-8859-1') as file2:
        file2.readline() 
        xml_soup = bs(file2.read(), 'xml')
        result = xml_soup.find_all('Disorder')
        orpha_sp={}
        for j in result:
            for i in j.find('OrphaNumber'):
                orpha= i.string
                result2 = j.find_all('ExternalReference')
                for i in result2:
                    for a in i.find('Source'):
                        if a.string == 'SwissProt':
                            for b in i.find('Reference'):
                                sp=b.string
                                orpha_sp[orpha]=orpha_sp.get(orpha,[])
                                orpha_sp[orpha].append(sp)
                                uniprot.append(sp)
    # print(len(uniprot))
    # print(len(set(uniprot)))
    return orpha_sp

def orphanet_poly(orph_swiss,humvar):
    poly = open(humvar)
    d_poly = {}
    orph_poly={}
    lpoly=[]
    ac=[]
    all_pol=[]
    for line in poly:
        if line[48:60] == 'Polymorphism':
            gene = line[0:9].strip()
            AC = line[10:20].strip()
            VAR = line[21:32].strip()
            AA_change=line[33:47].strip()
            dbSNP = line[62:76].strip()
            d_poly[AC] = d_poly.get(AC, {})
            d_poly[AC]['gene'] = d_poly[AC].get('gene', gene)
            d_poly[AC]['var'] = d_poly[AC].get('var', [])
            d_poly[AC]['var'].append(VAR)
            d_poly[AC]['snp'] = d_poly[AC].get('snp', [])
            d_poly[AC]['snp'].append(dbSNP)
            d_poly[AC]['aa_change'] = d_poly[AC].get('aa_change', [])
            d_poly[AC]['aa_change'].append(AA_change)

    for j in orph_swiss:
        for i in orph_swiss[j]:
           if i in d_poly:
               if [i,d_poly[i]['aa_change']] not in lpoly:
                    lpoly+=d_poly[i]['aa_change']
               orph_poly[j]=orph_poly.get(j, {})
               orph_poly[j]['AC']=orph_poly[j].get('AC',[])
               orph_poly[j]['AC'].append(i)
               ac.append(i)
               orph_poly[j]['var']=orph_poly[j].get('var',[])
               orph_poly[j]['var']+=d_poly[i]['var']
               
               orph_poly[j]['gene']=orph_poly[j].get('gene',d_poly[i]['gene'])

               orph_poly[j]['snp']=orph_poly[j].get('snp',[])
               orph_poly[j]['snp']+=d_poly[i]['snp']
               orph_poly[j]['aa_change']=orph_poly[j].get('aa_change',[])
               orph_poly[j]['aa_change']+=d_poly[i]['aa_change']
    
    print(len(set(lpoly)))
    return orph_poly,lpoly

def group_poly(orph_poly,group_orpha):
    g_poly={}
    all_poly=[]
    for i in group_orpha:
        for j in orph_poly:
            if j in group_orpha[i]:
                g_poly[i]= g_poly.get(i,{})
                g_poly[i]['AC']=g_poly[i].get('AC',[])
                g_poly[i]['AC']+=orph_poly[j]['AC']
                g_poly[i]['orpha']=g_poly[i].get('orpha',[])
                g_poly[i]['orpha'].append(j)
                g_poly[i]['var']=g_poly[i].get('var',[])
                g_poly[i]['var']+=orph_poly[j]['var']
                all_poly+=orph_poly[j]['var']
                g_poly[i]['snp']=g_poly[i].get('snp',[])
                g_poly[i]['snp']+=orph_poly[j]['snp']
                g_poly[i]['aa_change']=g_poly[i].get('aa_change',[])
                g_poly[i]['aa_change']+=orph_poly[j]['aa_change']
    for i in g_poly:
        print(i, len(set(g_poly[i]['var'])))
    return g_poly

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
    return l

def filter_cv(cv_cl):
    colnames = ['#AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinicalSignificance', 'ClinSigSimple', 'LastEvaluated', 'SNP', 'nsv/esv(dbVar)', 'RCVaccession', 'PhenotypeIDS', 'PhenotypeList','Origin','OrigininS',
    'Assembly','Chr1','Chr','Start','Stop','Ref','Alt',1,2,3,4,5,6,7,8]
    data = pd.read_csv(cv_cl, delimiter='\t', names=colnames,low_memory=False)
    d={}
    b={}
    c=0
    aorpha=[]
    aomim=[]
    avar=[]

    d_pred={}

    name = data.Name.tolist()
    ids = data.PhenotypeIDS.tolist()
    path=data.ClinicalSignificance.tolist()
    ch=data.Chr.tolist()
    start=data.Start.tolist()
    stop=data.Stop.tolist()
    ref=data.Ref.tolist()
    alt=data.Alt.tolist()
    gene=data.GeneSymbol.tolist()
    ass=data.Assembly.tolist()
    tp=data.Type.tolist()
    snp=data.SNP.tolist()

    for i in ids:
        c+=1
        n=name[c-1]
        if path[c-1]=='Pathogenic' and tp[c-1]=='single nucleotide variant' and ass[c-1]=='GRCh38':
            d[n]=d.get(n,{})
            d[n]['omim']=d[n].get('omim',[])
            d[n]['orphanet']=d[n].get('orphanet',[])
            d[n]['SNP']=d[n].get('SNP',[])
            d[n]['SNP'].append(snp[c-1])
            i=i.split(',')

            d_pred[n]=d_pred.get(n,{})
            d_pred[n]['path']=d_pred[n].get('path',path[c-1])
            d_pred[n]['gene']=d_pred[n].get('gene',gene[c-1])
            d_pred[n]['ch']=d_pred[n].get('ch',ch[c-1])
            d_pred[n]['start']=d_pred[n].get('start',start[c-1])
            d_pred[n]['stop']=d_pred[n].get('stop',stop[c-1])
            d_pred[n]['ref']=d_pred[n].get('ref',ref[c-1])
            d_pred[n]['alt']=d_pred[n].get('alt',alt[c-1])
            d_pred[n]['orphanet']=d_pred[n].get('orphanet',[])    


            for l in i:
                if l.find('|')!= -1:
                    i.remove(l)
                    i+=l.split('|')
            for j in i:
                o=j.find('OMIM') 
                if o != -1:
                    om=j[o+5: o+11]
                    d[n]['omim'].append(om)
                    aomim.append(om)
                r=j.find('Orphanet')
                if r != -1:
                    orph=j[r+14:]
                    d[n]['orphanet'].append(orph)
                    d_pred[n]['orphanet'].append(orph)
                    aorpha.append(orph)
                    
        
        if path[c-1]=='Benign' and tp[c-1]=='single nucleotide variant' and ass[c-1]=='GRCh38':
                b[n]=b.get(n,{})
                b[n]['omim']=b[n].get('omim',[])
                b[n]['orphanet']=b[n].get('orphanet',[])
                b[n]['SNP']=b[n].get('SNP',[])
                b[n]['SNP'].append(snp[c-1])
                avar.append(n)
                i=i.split(',')

                d_pred[n]=d_pred.get(n,{})
                d_pred[n]['path']=d_pred[n].get('path',path[c-1])
                d_pred[n]['gene']=d_pred[n].get('gene',gene[c-1])
                d_pred[n]['ch']=d_pred[n].get('ch',ch[c-1])
                d_pred[n]['start']=d_pred[n].get('start',start[c-1])
                d_pred[n]['stop']=d_pred[n].get('stop',stop[c-1])
                d_pred[n]['ref']=d_pred[n].get('ref',ref[c-1])
                d_pred[n]['alt']=d_pred[n].get('alt',alt[c-1])
                d_pred[n]['orphanet']=d_pred[n].get('orphanet',[])  

                for l in i:
                    if l.find('|')!= -1:
                        i.remove(l)
                        i+=l.split('|')
                for j in i:
                    o=j.find('OMIM') 
                    if o != -1:
                        om=j[o+5: o+11]
                        b[n]['omim'].append(om)
                        aomim.append(om)
                    r=j.find('Orphanet')
                    if r != -1:
                        orph=j[r+14:]
                        b[n]['orphanet'].append(orph)
                        d_pred[n]['orphanet'].append(orph)
                        aorpha.append(orph)  
    return d,b,d_pred
    
def group_clinvar(cv_filter,group_orpha):
    g_cv={}
    all_snp=[]
    snv=[]
    for i in group_orpha:
        g_cv[i]=g_cv.get(i,[])
        for j in cv_filter:
            if type(j)!= type(0.5):
                for n in cv_filter[j]['orphanet']:
                    if n in group_orpha[i]:
                        g_cv[i].append(j.replace(' ','_').replace('(','').replace(')',''))
                        # all_snp+=cv_filter[j]['SNP']
    for i in g_cv:
        print(i,len(set(g_cv[i])))
    print(len(g_cv))
    # print(g_cv)
    return g_cv





if __name__ == '__main__':
    humvar = sys.argv[1]
    xml1 = sys.argv[2]
    group_csv = sys.argv[3]
    xml2=sys.argv[4]
    cv_cl=sys.argv[5]
    #mim_var,lvar,gene_up = omim_var(humvar)
    # orphanet_mim = orphanet_omim(xml1)
    orph, xml = orphanet_group(group_csv)
    group_orpha = group_orphanet(orph, xml)
    # orpha_var, all_var= orphanet_var(mim_var, orphanet_mim)
    #group_var=group_var(orpha_var, group_orpha)
    # orph_swiss=orphanet_swiss(xml2)
    # orph_poly,poly=orphanet_poly(orph_swiss,humvar)
    # group_polym = group_poly(orph_poly,group_orpha)
    cv_filter_p,cv_filter_b,pred= filter_cv(cv_cl)
    group_cv = group_clinvar(cv_filter_b,group_orpha)
    #group_cv_pred= group_clinvar(pred,group_orpha)

    # with open('omim_var.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in list(mim_var.keys()):
    #         c = ','
    #         mim_var[i]['var'] = c.join(mim_var[i]['var'])
    #         mim_var[i]['snp'] = c.join(mim_var[i]['snp'])
    #         all_writer.writerow(
    #             [i, mim_var[i]['omim'], mim_var[i]['gene'], mim_var[i]['var'], mim_var[i]['snp']])

    # with open('orphanet_mim.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in list(orphanet_mim.keys()):
    #         c = ','
    #         all_writer.writerow(
    #             [i, c.join(orphanet_mim[i][0]), orphanet_mim[i][1]])

    # with open('group_orphanet.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in list(group_orpha.keys()):
    #         c = ','
    #         all_writer.writerow([int(i), c.join(set(group_orpha[i]))])

    # with open('orphanet_var2.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in list(orpha_var.keys()):
    #         c = ','
    #         all_writer.writerow([i,c.join(set(orpha_var[i]['gene'])),c.join(set(orpha_var[i]['change']))])
 
    # with open('group_var.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in list(group_var.keys()):
    #         c = ','
    #         all_writer.writerow([int(i), c.join(set(group_var[i]['orpha'])),c.join(set(group_var[i]['val'])),c.join(set(group_var[i]['var'])),c.join(set(group_var[i]['snp']))])

    # with open('orphanet_swiss.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in orph_swiss:
    #         c = ','
    #         all_writer.writerow([i, c.join(set(orph_swiss[i]))])

    # with open('orphanet_poly.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in list(orph_poly.keys()):
    #         c = ','
    #         all_writer.writerow([i, c.join(set(orph_poly[i]['AC'])),c.join(set(orph_poly[i]['var'])),c.join(set(orph_poly[i]['snp'])),c.join(set(orph_poly[i]['aa_change']))])

    # with open('group_poly.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in list(group_polym.keys()):
    #         c = ','
    #         all_writer.writerow([int(i),c.join(set(group_polym[i]['orpha'])),c.join(set(group_polym[i]['AC'])),c.join(set(group_polym[i]['var'])),c.join(set(group_polym[i]['snp'])),c.join(set(group_polym[i]['aa_change']))])

    # with open('cv_snv_patho.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     c=','
    #     for i in cv_filter:
    #         all_writer.writerow([i,c.join(set(cv_filter[i]['omim'])),c.join(set(cv_filter[i]['orphanet'])),c.join(set(cv_filter[i]['SNP']))])

    # with open('cv_snv_benign.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     c=','
    #     cv_filter=cv_filter_b
    #     for i in cv_filter:
    #         all_writer.writerow([i,c.join(set(cv_filter[i]['omim'])),c.join(set(cv_filter[i]['orphanet'])),c.join(set(cv_filter[i]['SNP']))])

    # with open('cv_snv_patho_benign.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     cv_filter=pred
    #     for i in cv_filter:
    #         all_writer.writerow([i,cv_filter[i]['path'],cv_filter[i]['gene'],cv_filter[i]['ch'],cv_filter[i]['start'],cv_filter[i]['ref'],cv_filter[i]['alt']])


    # with open('group_clinvar.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in group_cv.keys():
    #         c = ','
    #         all_writer.writerow([int(i), c.join(set(group_cv[i]))])

    # with open('group_clinvar_pred.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in group_cv_pred.keys():
    #         c = ','
    #         all_writer.writerow([i, c.join(set(group_cv_pred[i]))])
    #         #all_writer.writerow([i, set(group_cv_pred[i])])
    

    # with open('omimvar.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     for i in lvar:
    #         all_writer.writerow([i[0],i[1]])


    # with open('orphavar2.csv', mode='w') as all_file:
    #     all_writer = csv.writer(all_file)
    #     c=[]
    #     for i in range(len(all_var)-1):
    #         for j in all_var[i][0]:
    #             if [all_var[i][1],j,all_var[i][2]] not in c: 
    #                 c.append([all_var[i][1],j,all_var[i][2]])
    #     for l in c: 
    #         all_writer.writerow([l[1]+'-'+l[0],l[2]])


    # with open('orphapoly2.csv', mode='w') as all_file: #orph_poly
    #     all_writer = csv.writer(all_file)
    #     c=[]
    #     for i in orph_poly:
    #         for j in orph_poly[i]['aa_change']:
    #             if [j+'-'+orph_poly[i]['gene'],i] not in c:
    #                 c.append([j+'-'+orph_poly[i]['gene'],i])
    #                 all_writer.writerow([j+'-'+orph_poly[i]['gene'],i])

    

