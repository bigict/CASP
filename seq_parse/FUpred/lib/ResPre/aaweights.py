# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 22:39:40 2017

@author: 74297
"""

import numpy as np
import os
from subprocess import Popen, PIPE, STDOUT
from numba import jit
from io import BytesIO
aadic={
        'A':1,
        'B':0,
        'C':2,
        'D':3,
        'E':4,
        
        'F':5,
        'G':6,
        'H':7,
        'I':8,
        'J':0,
        
        'K':9,
        'L':10,
        'M':11,
        'N':12,
        'O':0, 
        
        'P':13,
        'Q':14,
        'R':15,
        'S':16,
        'T':17,
        'U':0,
        
        'V':18,
        'W':19,
        'X':0,
        'Y':20,
        'Z':0,
        '-':0,
        '*':0,
        }


def getb62(bfile,aawightsof_):
    
    m=np.zeros([21,21])
    lines=open(bfile).readlines()
    aaindexs=lines[0].split()
    for i in range(1,len(lines)):
        segs=lines[i].split()
        aa=segs[0]
        for j in range(1,len(segs)):
            m[aadic[aaindexs[j-1]]][aadic[aa]]=float(segs[j])
    for i in range(21):#set - to 0 which is no effect
        m[0,i]=m[i,0]=aawightsof_
    return m
def readsequence_old(seq_file):
    lines=open(seq_file).readlines()
    seq=lines[1].strip()
    aas=[aadic[aa] for aa in seq]
    return aas
def readsequence(seq_file):
    lines=open(seq_file).readlines()
    lines=[line.strip() for line in lines]
    seq=lines[1].strip()  
    aalines=''
    for i in range(1,len(lines)):
        aalines+=lines[i]
    aas=[aadic[aa] for aa in aalines]
    return aas
def getweights_(msafile,seq_id):
    exefile='bin/colors_weight'
    cmd=exefile+' '+msafile+' '+str(seq_id)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    return np.genfromtxt(p.stdout)
def getweights(msafile,seq_id):
    exefile='bin/colors_weight_win.exe'
    cmd=exefile+' '+msafile+' '+str(seq_id)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    output,error = p.communicate()
    #output=output.decode("utf-8")
    return np.genfromtxt(BytesIO(output)).flatten()
def getweights_out(msafile,seq_id,outfile):
    exefile='colors_weight_win.exe'
    cmd=exefile+' '+msafile+' '+str(seq_id)+' >'+outfile
    os.system(cmd)
    #p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    #output,error = p.communicate()
    #print(output)
    #return np.genfromtxt(output)    
def henikoff_weight(msa):
    n,p=msa.shape
    unique_sum=np.zeros(p)
    for i in range(p):
        myset=set(msa[:,i])
        #myset.discard(0)
        unique_sum[i]=len(list(myset))
    count_dict=np.zeros([21,p])
    for i in range(21):
        for j in range(p):
            count_dict[i,j]= list(msa[:,j]).count(i)
    weights=np.zeros(n)
    for i in range(n):
        aweight=0.0
        for j in range(p):
            #if msa[i,j]!=0:
           aweight+=1.0/(unique_sum[j]*count_dict[msa[i,j],j])
        weights[i]=aweight
    return weights
            
def read_msa(file_path):    
    lines=open(file_path).readlines()  
    lines=[line.strip() for line in lines]   
    n=len(lines)
    d=len(lines[0]) #CR AND LF  
    msa=np.zeros([n,d],dtype=int)
    for i in  range(n):
        aline=lines[i]
        for j in range(d):
            msa[i,j]=aadic[aline[j]]
    return msa
def generateweights(seq,aln,bm,weight1,alpha=0.0):
    #weights=exp(alpha*bm)
    [n,p]=aln.shape
    bm=np.exp(alpha*bm)  
    weights=np.zeros([n,p]);
    
    for i in range(n):
        for j in range(p):
            #print(weight1[i],seq[j],aln[i,j])
            weights[i,j]=weight1[i]/bm[seq[j],aln[i,j]]
    return weights
def cal_large_matrix_newweights(msa,weight):
    ALPHA=21
    pseudoc=1
    M=msa.shape[0]
    N=msa.shape[1]
    pab=np.zeros((ALPHA,ALPHA))
    pa=np.zeros((N,ALPHA))
    cov=np.zeros((N*ALPHA,N*ALPHA))
    for i in range(N):
        for aa in range(ALPHA):
            pa[i,aa] = pseudoc
        neff=0.0
        for k in range(M):
            pa[i,msa[k,i]]+=weight[k,i]
            neff+=weight[k,i]
        for aa in range(ALPHA):
            pa[i,aa] /=pseudoc * ALPHA * 1.0 + neff
    #print(pab)
    for i in range(N):
        for j in range(i,N):
            for a in range(ALPHA):
                for b in range(ALPHA):
                    if i ==j :
                        if a==b :
                            pab[a,b]=pa[i,a]
                        else:
                            pab[a,b]=0.0
                    else:
                        pab[a,b] = pseudoc *1.0 /ALPHA
            if(i!=j):
                neff2=0;
                for k in range(M):
                    a=msa[k,i]
                    b=msa[k,j]
                    tmp=np.sqrt(weight[k,i]*weight[k,j])
                    pab[a,b]+=tmp
                    neff2+=tmp
                for a in range(ALPHA):
                    for b in range(ALPHA):
                        pab[a,b] /= pseudoc*ALPHA*1.0 +neff2
            for a in range(ALPHA):
                for b in range(ALPHA):
                    if(i!=j or a==b):
                        if (pab[a][b] > 0.0):
                            cov[i*21+a][j*21+b]=pab[a][b] - pa[i][a] * pa[j][b]
                            cov[j*21+b][i*21+a]=cov[i*21+a][j*21+b]
            #cov1=cal_cov(pab,pa[<unsigned int>i],pa[<unsigned int>j])
            #cov[<unsigned int>i,<unsigned int>j]=cov[<unsigned int>j,<unsigned int>i]=cov1
    #print(pab)
    #print(pa)
    return cov
@jit   
def cal_large_matrix1(msa,weight):
    #output:21*l*21*l
    ALPHA=21
    pseudoc=1
    M=msa.shape[0]
    N=msa.shape[1]
    pab=np.zeros((ALPHA,ALPHA))
    pa=np.zeros((N,ALPHA))
    cov=np.zeros([N*ALPHA,N*ALPHA ])
    for i in range(N):
        for aa in range(ALPHA):
            pa[i,aa] = pseudoc
        neff=0.0
        for k in range(M):
            pa[i,msa[k,i]]+=weight[k]
            neff+=weight[k]
        for aa in range(ALPHA):
            pa[i,aa] /=pseudoc * ALPHA * 1.0 + neff
    #print(pab)
    for i in range(N):
        for j in range(i,N):
            for a in range(ALPHA):
                for b in range(ALPHA):
                    if i ==j :
                        if a==b :
                            pab[a,b]=pa[i,a]
                        else:
                            pab[a,b]=0.0
                    else:
                        pab[a,b] = pseudoc *1.0 /ALPHA
            if(i!=j):
                neff2=0;
                for k in range(M):
                    a=msa[k,i]
                    b=msa[k,j]
                    tmp=weight[k]
                    pab[a,b]+=tmp
                    neff2+=tmp
                for a in range(ALPHA):
                    for b in range(ALPHA):
                        pab[a,b] /= pseudoc*ALPHA*1.0 +neff2
            for a in range(ALPHA):
                for b in range(ALPHA):
                    if(i!=j or a==b):
                        if (pab[a][b] > 0.0):
                            cov[i*21+a][j*21+b]=pab[a][b] - pa[i][a] * pa[j][b]
                            cov[j*21+b][i*21+a]=cov[i*21+a][j*21+b]
                            #cov[a*ALPHA+b,i,j,]=pab[a][b] - pa[i][a] * pa[j][b]
                            #cov[a*ALPHA+b,j,i,]=cov[a*ALPHA+b,i,j]
            #cov1=cal_cov(pab,pa[<unsigned int>i],pa[<unsigned int>j])
            #cov[<unsigned int>i,<unsigned int>j]=cov[<unsigned int>j,<unsigned int>i]=cov1
    #print(pab)
    #print(pa)
    return cov 
@jit
def cal_large_matrix2(msa,weight):
    #output:441*l*l
    ALPHA=21
    pseudoc=1
    M=msa.shape[0]
    N=msa.shape[1]
    pab=np.zeros((ALPHA,ALPHA))
    pa=np.zeros((N,ALPHA))
    cov=np.zeros([ALPHA*ALPHA,N,N])
    for i in range(N):
        for aa in range(ALPHA):
            pa[i,aa] = pseudoc
        neff=0.0
        for k in range(M):
            pa[i,msa[k,i]]+=weight[k]
            neff+=weight[k]
        for aa in range(ALPHA):
            pa[i,aa] /=pseudoc * ALPHA * 1.0 + neff
    #print(pab)
    for i in range(N):
        for j in range(i,N):
            for a in range(ALPHA):
                for b in range(ALPHA):
                    if i ==j :
                        if a==b :
                            pab[a,b]=pa[i,a]
                        else:
                            pab[a,b]=0.0
                    else:
                        pab[a,b] = pseudoc *1.0 /ALPHA
            if(i!=j):
                neff2=0;
                for k in range(M):
                    a=msa[k,i]
                    b=msa[k,j]
                    tmp=weight[k]
                    pab[a,b]+=tmp
                    neff2+=tmp
                for a in range(ALPHA):
                    for b in range(ALPHA):
                        pab[a,b] /= pseudoc*ALPHA*1.0 +neff2
            for a in range(ALPHA):
                for b in range(ALPHA):
                    if(i!=j or a==b):
                        if (pab[a][b] > 0.0):
                            #cov[i*21+a][j*21+b]=pab[a][b] - pa[i][a] * pa[j][b]
                            #cov[j*21+b][i*21+a]=cov[i*21+a][j*21+b]
                            cov[a*ALPHA+b,i,j,]=pab[a][b] - pa[i][a] * pa[j][b]
                            cov[a*ALPHA+b,j,i,]=cov[a*ALPHA+b,i,j]
            #cov1=cal_cov(pab,pa[<unsigned int>i],pa[<unsigned int>j])
            #cov[<unsigned int>i,<unsigned int>j]=cov[<unsigned int>j,<unsigned int>i]=cov1
    #print(pab)
    #print(pa)
    return cov    

if __name__ == "__main__":
    import mycex
    import scipy.misc
    id='psicov/aln/1hxrA.aln'
    msa=mycex.read_msa(id)
    weights=henikoff_weight(msa)
    print(weights)
    