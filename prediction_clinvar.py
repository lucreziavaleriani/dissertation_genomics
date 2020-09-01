#!usr/bin/python3
import sys
import shutil
import tempfile
import urllib.request
import ssl
import csv
import pandas as pd
from sklearn import metrics
import requests
from bs4 import  BeautifulSoup as bs
ssl._create_default_https_context = ssl._create_unverified_context

aa={'K':'Lys','R':'Arg','H':'His','D':'Asp','E':'Glu','A':'Ala','V':'Val','I':'Ile','L':'Leu',
'M':'Met','F':'Phe','Y':'Tyr','W':'Trp','C':'Cys','U':'Sec','G':'Gly','P':'Pro','S':'Ser','T':'Thr',
'N':'Asn','Q':'Gln'}

#benign = poly = 1
#patho = disease = 0

def file_creation(clinvar):
    colnames = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'EFFECT','GENE','MUT','VAR', 'CODING', 'PREDICTION', 'SCORE', 'FDR', 'PhyloP100', 'AvgPhyloP100']
    all_d={}
    data = pd.read_csv(clinvar, delimiter='\t', names=colnames)
       
    ch=data.CHROM.tolist()
    pos=data.POS.tolist()
    L_id=data.ID.tolist()
    ref=data.REF.tolist()
    alt=data.ALT.tolist()
    expected=data.EFFECT.tolist()
    gene=data.GENE.tolist()
    mut=data.MUT.tolist()
    var=data.VAR.tolist()
    code=data.CODING.tolist()
    score=data.SCORE.tolist() #if score=0.5 --> benign 
    pred=data.PREDICTION.tolist()

    c=0
    for i in var[1:]:
        c+=1
        all_d[i]=all_d.get(i,{})
        all_d[i]['pred']=pred[c]
        all_d[i]['exp']=expected[c]
        all_d[i]['score']=score[c]

    exp=[]
    for i in all_d:
        exp.append(all_d[i]['exp'])

    return all_d,exp

#now analyze the prediction in each group
#riassegnare ad ogni variante il proprio orphanet


def group(group,all_d):
    colnames=['group','cv']
    data = pd.read_csv(group, delimiter=',', names=colnames)
    
    group=data.group.tolist()
    cv=data.cv.tolist()

    group_var={}
    c=0
    for n in cv:
        if type(n) != type(0.5):
            i=n.split(',')
            c+=1            
            for k in all_d:
                if type(k) != type(0.5):
                    if k.replace('(','').replace(')','') in i:
                        group_var[group[c-1]]=group_var.get(group[c-1],[])
                        group_var[group[c-1]].append(k)
    return group_var


def group_stats(group_var,d):
    bs={}
    all_exp={}
    for i in group_var:
        all_exp[i]=all_exp.get(i,[])
        for j in group_var[i]:
                all_exp[i].append(d[j]['exp'])

    all_thr={}
    for i in group_var:
        all_thr[i] = threshold(group_var[i],d)
    
    for i in all_exp:
        bs[i]=result(all_exp[i],all_thr[i])
        #best_stat[i]=best_stats(bs)
    
    print(bs.keys())
    return bs


def analysis(pred,exp):
    cm=metrics.confusion_matrix(exp, pred,labels=['Pathogenic','Benign'])
    acc=metrics.accuracy_score(exp,pred)
    mcc= metrics.matthews_corrcoef(exp,pred)
    #auc= metrics.roc_auc_score(exp,pred)
    tn, fp, fn, tp = metrics.confusion_matrix(exp, pred,labels=['Pathogenic','Benign']).ravel()
    ppv=tp/(tp+fp)
    npv=tn/(fn+tn)
    d= {'cm':cm,'acc':acc,'mcc':mcc,'tpr':tp/(tp+fn),'tnr':tn/(tn+fp),'ppv':ppv,'npv':npv}
    return d

def threshold(new,d):
    new_pred={}
    thr=[0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9]
    for i in thr:
        pred=[]
        for j in new:
            if float(d[j]['score']) > i:
                pred.append('Pathogenic')
            if float(d[j]['score']) <= i:
                pred.append('Benign')
        new_pred[i]=pred
    return new_pred

def result(exp,new_pred):
    metrics={}
    for i in new_pred:
        metrics[i]=analysis(new_pred[i],exp)
    for i in metrics:
        print(i,round(metrics[i]['acc'],5),round(metrics[i]['ppv'],5),'&',round(metrics[i]['npv'],5))
    return metrics


if __name__ == '__main__':
    clinvar = sys.argv[1]
    group_cv = sys.argv[2]

    d,exp = file_creation(clinvar)
    # th = threshold(d)
    # res = result(exp,th)
    g_var=group(group_cv,d)
    g_stats = group_stats(g_var,d)



# with open('clinvar_pred_group.csv', mode='w') as all_file:
#         all_writer = csv.writer(all_file)
#         cv_filter = g_stats
#         for i in cv_filter:
#             for j in cv_filter[i]:
#                 all_writer.writerow([i,j,cv_filter[i][j]['acc'],cv_filter[i][j]['mcc'],cv_filter[i][j]['tpr'],cv_filter[i][j]['tnr'],cv_filter[i][j]['ppv'],cv_filter[i][j]['npv']])
