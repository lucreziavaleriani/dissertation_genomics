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

def file_creation(pred_p,pred_d):
    colnames = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'PROT', 'CODING', 'PREDICTION', 'SCORE', 'FDR', 'PhyloP100', 'AvgPhyloP100']
    c=[pred_p,pred_d]
    new={}
    new_g={}
    all_pred=[]
    all_exp=[]
    for i in c:
        data = pd.read_csv(i, delimiter='\t', names=colnames)
       
        ch=data.CHROM.tolist()
        pos=data.POS.tolist()
        l_id=data.ID.tolist()
        ref=data.REF.tolist()
        alt=data.ALT.tolist()
        prot=data.PROT.tolist()
        code=data.CODING.tolist()
        score=data.SCORE.tolist() #if score=0.5 --> benign 
        prediction=data.PREDICTION.tolist()
        
        if i == 'poly.csv':
            exp= ['Benign']*(len(prediction))
        else:
            exp= ['Pathogenic']*(len(prediction))
        
        for n in range(1,len(prediction)):
            new[prot[n]+'-'+pos[n]+'-'+ch[n]]=new.get(pos[n]+'-'+pos[n]+'-'+ch[n],{})
            new[prot[n]+'-'+pos[n]+'-'+ch[n]]['pred']=prediction[n]
            new[prot[n]+'-'+pos[n]+'-'+ch[n]]['exp']=exp[n]
            new[prot[n]+'-'+pos[n]+'-'+ch[n]]['score']=score[n]

            p=prot[n].find(':')
            if p!=-1:
                gene=prot[n][:p]
                var=prot[n][p+1:]
                s=aa[var[0]]
                e=aa[var[-1]]
                new_var='p.'+s+var[1:-1]+e
                new_g[new_var+'-'+gene]=new_g.get(new_var+'-'+gene,{})
                new_g[new_var+'-'+gene]['pred']=prediction[n]
                new_g[new_var+'-'+gene]['exp']=exp[n]
                new_g[new_var+'-'+gene]['score']=score[n]
    
    for i in new:
        all_exp.append(new[i]['exp'])
        all_pred.append(new[i]['pred'])
    return new,all_pred,all_exp,new_g

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

def threshold(new):
    new_pred={}
    thr=[0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9]
    for i in thr:
        pred=[]
        for j in new:
            if float(new[j]['score']) > i:
                pred.append('Pathogenic')
            if float(new[j]['score']) <= i:
                pred.append('Benign')
        new_pred[i]=pred
    return new_pred


def result(exp,new_pred):
    metrics={}
    for i in new_pred:
        metrics[i]=analysis(new_pred[i],exp)
    for i in metrics:
        print(i,round(metrics[i]['acc'],5),'&',round(metrics[i]['mcc'],5),'&',round(metrics[i]['tpr'],5),'&',round(metrics[i]['tnr'],5),'&',round(metrics[i]['ppv'],5),'&',round(metrics[i]['npv'],5),'\\')
    return metrics

def best_stats(metrics):
    best=''
    for i in metrics:
        if best=='':
            best=metrics[i]
        elif metrics[i]['acc'] > best['acc']:
            best= metrics[i]['acc']
    return best

#now analyze the prediction in each group
#riassegnare ad ogni variante il proprio orphanet

def orphavar(g,orpha_var):
    colnames=['NAME','ORPHA']
    data = pd.read_csv(orpha_var, delimiter=',', names=colnames)
    name=data.NAME.tolist()
    orpha=data.ORPHA.tolist()

    all_name=[]

    o_v={}
    for i in g:
        c=0
        for j in name:
            c+=1
            if i == j:
                o_v[orpha[c-1]]=o_v.get(orpha[c-1],[])
                o_v[orpha[c-1]].append(j)
    return o_v

def group(group,o_v):
    colnames=['group','orpha']
    if group== 'group_poly2.csv':
        data = pd.read_csv(group, delimiter=';', names=colnames)
    else: 
        data = pd.read_csv(group, delimiter=',', names=colnames)
    
    group=data.group.tolist()
    orpha=data.orpha.tolist()

    group_var={}
    c=0
    for n in orpha:
        i=n.split(',')
        c+=1            
        for j in i:
            for k in o_v:
                if int(k)==int(j):
                    group_var[group[c-1]]=group_var.get(group[c-1],[])
                    group_var[group[c-1]]+=o_v[k]
    return group_var


def group_stats(group_var,group_poly,g):
    bs={}
    
    d={}

    all_exp={}
    for i in group_var: # all_exp[i]=all_exp.get(i,[])
        d[i]=d.get(i,{})
        for j in group_var[i]:#all_exp[i].append(g[j]['exp'])
            d[i][j]=g[j]

    for i in group_poly:
        d[i]=d.get(i,{})
        for j in group_poly[i]:
            d[i][j]=g[j]
    
    for j in d:
        all_exp[j]=all_exp.get(i,[])
        for i in d[j]:
            all_exp[j].append(g[i]['exp'])

    all_thr={}
    for i in d:
        all_thr[i] = threshold(d[i])
    
    for i in all_exp:
        bs[i]=result(all_exp[i],all_thr[i])
        #best_stat[i]=best_stats(bs)
    
    #print(best_stat)
    return bs


if __name__ == '__main__':
    disease = sys.argv[1]
    poly = sys.argv[2]
    group_var = sys.argv[3]
    orpha_var=sys.argv[4]
    orpha_poly=sys.argv[5]
    group_poly=sys.argv[6]

    new,pred,exp,g= file_creation(poly,disease)
    
    orp_v=orphavar(g,orpha_var)
    orp_p=orphavar(g,orpha_poly)
    
    group_v=group(group_var,orp_v)
    group_p=group(group_poly,orp_p)


    stat=group_stats(group_v,group_p,g)
    
    # thr=threshold(new)
    # stat2=result(exp,thr)

# with open('analysis_group_pred.csv', mode='w') as all_file:
#         all_writer = csv.writer(all_file)
#         cv_filter = stat
#         for i in cv_filter:
#             for j in cv_filter[i]:
#                 all_writer.writerow([i,j,cv_filter[i][j]['acc'],cv_filter[i][j]['mcc'],cv_filter[i][j]['tpr'],cv_filter[i][j]['tnr'],cv_filter[i][j]['ppv'],cv_filter[i][j]['npv']])
