#!/usr/bin/python
from __future__ import print_function
import getopt
import sys
import os
import numpy as np
import pandas as pd
import math
import itertools
from itertools import repeat
import argparse
import csv
from collections import Counter
import re
import glob
import time
from time import sleep
from tqdm import tqdm
from argparse import RawTextHelpFormatter
import uuid
import warnings
warnings.filterwarnings("ignore") 
std = list("ACDEFGHIKLMNPQRSTVWY")
def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
paths = [find('aa_attr_group.csv','/')]
paths_1 = paths[0].split('/')[:-1]
paths_2 = '/'.join(paths_1)
PCP= pd.read_csv(paths_2+'/PhysicoChemical.csv', header=None)
AAindices = paths_2+'/aaind.txt'
AAIndex = pd.read_csv(paths_2+'/aaindex.csv',index_col='INDEX');
AAIndexNames = pd.read_csv(paths_2+'/AAIndexNames.csv',header=None);
def aac_comp(file,out):
    f = open(out, 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    zz = df.iloc[:,0]
    print("AAC_A,AAC_C,AAC_D,AAC_E,AAC_F,AAC_G,AAC_H,AAC_I,AAC_K,AAC_L,AAC_M,AAC_N,AAC_P,AAC_Q,AAC_R,AAC_S,AAC_T,AAC_V,AAC_W,AAC_Y,")
    for j in zz:
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            print("%.2f"%composition, end = ",")
        print("")
    f.truncate()
def aac_wp(file,out):
     readseq(file,'input_sam.csv')
     aac_comp('input_sam.csv','tempfile_out')
     df = pd.read_csv('tempfile_out')
     df.iloc[:,:-1].to_csv(out,index=None)
     os.remove('tempfile_out')
     os.remove('input_sam.csv')
def aac_nt(file,out,n):
     readseq(file,'input_sam.csv')
     file1 = nt('input_sam.csv',n)
     file1.to_csv('sam_input.csv', index=None, header=False)
     aac_comp('sam_input.csv','tempfile_out')
     df = pd.read_csv('tempfile_out')
     df.columns = 'N'+df.columns
     df.iloc[:,:-1].to_csv(out,index=None)
     os.remove('sam_input.csv')
     os.remove('tempfile_out')
     os.remove('input_sam.csv')
def aac_ct(file,out,c):
     readseq(file,'input_sam.csv')
     file1 = ct('input_sam.csv',c)
     file1.to_csv('sam_input.csv', index=None, header=False)
     aac_comp('sam_input.csv','tempfile_out')
     df = pd.read_csv('tempfile_out')
     df.columns = 'C'+df.columns
     df.iloc[:,:-1].to_csv(out,index=None)
     os.remove('sam_input.csv')
     os.remove('tempfile_out')
     os.remove('input_sam.csv')
def aac_rt(file,out,n,c):
    readseq(file,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aac_comp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.iloc[:,:-1].to_csv(out,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def aac_nct(file,out,n):
    readseq(file,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aac_comp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.iloc[:,:-1].to_csv(out,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def aac_st(file,out,v):
    readseq(file,'input_sam.csv')
    file1 = split('input_sam.csv',v)
    file2 = file1.iloc[:,0]
    file2.to_csv('sam_input.csv', index=None, header=None)
    aac_comp('sam_input.csv','tempfile_out')
    df4_1 = pd.read_csv("tempfile_out")
    df4 = df4_1.iloc[:,:-1]
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*20)):
        pp.append(ss[i:i+(v*20)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
def dpc_comp(file,q,out):
    filename, file_extension = os.path.splitext(file)
    f = open(out, 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    zz = df1.iloc[:,0]
    for s in std:
        for u in std:
            print("DPC"+str(q)+"_"+s+u, end=',')
    print("")
    for i in range(0,len(zz)):
        for j in std:
            for k in std:
                count = 0
                temp = j+k
                for m3 in range(0,len(zz[i])-q):
                    b = zz[i][m3:m3+q+1:q]
                   # b.upper()

                    if b == temp:
                        count += 1
                    composition = (count/(len(zz[i])-(q)))*100
                print("%.2f" %composition, end = ',')
        print("")
    f.truncate()
def dpc_wp(seq,result_filename,lg):
    readseq(seq,'input_sam.csv')
    dpc_comp('input_sam.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def dpc_nt(seq,result_filename,n,lg):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    dpc_comp('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def dpc_ct(seq,result_filename,c,lg):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    dpc_comp('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def dpc_rt(seq,result_filename,n,c,lg):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    dpc_comp('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def dpc_nct(seq,result_filename,n,lg):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    dpc_comp('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def dpc_st(seq,result_filename,sp,lg):
    readseq(seq,'input_sam.csv')
    dpc_split('input_sam.csv',lg,sp,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def dpc_split(file,q,n,out):
    filename,file_ext = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    k1 = []
    for w in range(0,len(df2)):
        s = 0
        k2 = []
        r = 0
        if len(df2[0][w])%n == 0:
            k2.extend(repeat(int(len(df2[0][w])/n),n))
        else:
            r = int(len(df2[0][w])%n)
            k2.extend(repeat(int(len(df2[0][w])/n),n-1))
            k2.append((int(len(df2[0][w])/n)+r))
        for j in k2:
            df3 = df2[0][w][s:j+s]
            k1.append(df3)
            s = j+s
    f = open(out, 'w')
    sys.stdout = f
    for h in range(1,n+1):
        for e in std:
            for r in std:
                print('DPC'+str(q)+'_'+e+r+'_s'+str(h), end=",")
    print("")
    for i in range(0,len(k1),n):
        k4 = k1[i:i+n]
        for j in k4:
            for k in std:
                for l in std:
                    count = 0
                    temp = k+l
                    for m3 in range(0,(len(j)-q)):
                        b = j[m3:m3+q+1:q]
                        if b == temp:
                            count += 1
                        composition = (count/(len(j)-q))*100
                    print("%.2f" %composition, end = ',')
        print("")
    f.truncate()

#############################TPC_COMP##############################
def tpc_comp(file,out):
    f = open(out, 'w')
    sys.stdout = f
    df = pd.read_csv(file, header = None)
    zz = df.iloc[:,0]
    for s in std:
       for u in std:
           for m in std:
               print('TPC_'+s+u+m, end=",")
    print("")
    for i in range(0,len(zz)):
        for j in std:
            for k in std:
                for m1 in std:
                    count = 0
                    temp = j+k+m1
                    for m3 in range(0,len(zz[i])):
                        b = zz[i][m3:m3+3]
                        if b == temp:
                            count += 1
                        composition = (count/(len(zz[i])-2))*100
                    print("%.2f" %composition, end = ',')
        print("")
    f.truncate()
##################################################################################
def tpc_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    tpc_comp('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def tpc_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    tpc_comp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def tpc_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    tpc_comp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def tpc_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    tpc_comp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def tpc_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    tpc_comp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def tpc_st(file,out,v):
    readseq(file,'input_sam.csv')
    file1 = split('input_sam.csv',v)
    file2 = file1.iloc[:,0]
    file2.to_csv('sam_input.csv', index=None, header=None)
    tpc_comp('sam_input.csv','tempfile_out')
    df4_1 = pd.read_csv("tempfile_out")
    df4 = df4_1.iloc[:,:-1]
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*8000)):
        pp.append(ss[i:i+(v*8000)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
###########################atom###############
def atc(file,out):
    atom=pd.read_csv(paths_2+"/atom.csv",header=None)
    at=pd.DataFrame()
    i = 0
    C_atom = []
    H_atom = []
    N_atom = []
    O_atom = []
    S_atom = []

    while i < len(atom):
        C_atom.append(atom.iloc[i,1].count("C"))
        H_atom.append(atom.iloc[i,1].count("H"))
        N_atom.append(atom.iloc[i,1].count("N"))
        O_atom.append(atom.iloc[i,1].count("O"))
        S_atom.append(atom.iloc[i,1].count("S"))
        i += 1
    atom["C_atom"]=C_atom
    atom["O_atom"]=O_atom
    atom["H_atom"]=H_atom
    atom["N_atom"]=N_atom
    atom["S_atom"]=S_atom
##############read file ##########
    test1 = pd.read_csv(file,header=None)
    dd = []
    for i in range(0, len(test1)):
        dd.append(test1[0][i].upper())
    test = pd.DataFrame(dd)
    count_C = 0
    count_H = 0
    count_N = 0
    count_O = 0
    count_S = 0
    count = 0
    i1 = 0
    j = 0
    k = 0
    C_ct = []
    H_ct = []
    N_ct = []
    O_ct = []
    S_ct = []
    while i1 < len(test) :
        while j < len(test[0][i1]) :
            while k < len(atom) :
                if test.iloc[i1,0][j]==atom.iloc[k,0].replace(" ","") :
                    count_C = count_C + atom.iloc[k,2]
                    count_H = count_H + atom.iloc[k,3]
                    count_N = count_N + atom.iloc[k,4]
                    count_O = count_O + atom.iloc[k,5]
                    count_S = count_S + atom.iloc[k,6]
                #count = count_C + count_H + count_S + count_N + count_O
                k += 1
            k = 0
            j += 1
        C_ct.append(count_C)
        H_ct.append(count_H)
        N_ct.append(count_N)
        O_ct.append(count_O)
        S_ct.append(count_S)
        count_C = 0
        count_H = 0
        count_N = 0
        count_O = 0
        count_S = 0
        j = 0
        i1 += 1
    test["C_count"]=C_ct
    test["H_count"]=H_ct
    test["N_count"]=N_ct
    test["O_count"]=O_ct
    test["S_count"]=S_ct

    ct_total = []
    m = 0
    while m < len(test) :
        ct_total.append(test.iloc[m,1] + test.iloc[m,2] + test.iloc[m,3] + test.iloc[m,4] + test.iloc[m,5])
        m += 1
    test["count"]=ct_total
##########final output#####
    final = pd.DataFrame()
    n = 0
    p = 0
    C_p = []
    H_p = []
    N_p = []
    O_p = []
    S_p = []
    while n < len(test):
        C_p.append((test.iloc[n,1]/test.iloc[n,6])*100)
        H_p.append((test.iloc[n,2]/test.iloc[n,6])*100)
        N_p.append((test.iloc[n,3]/test.iloc[n,6])*100)
        O_p.append((test.iloc[n,4]/test.iloc[n,6])*100)
        S_p.append((test.iloc[n,5]/test.iloc[n,6])*100)
        n += 1
    final["ATC_C"] = C_p
    final["ATC_H"] = H_p
    final["ATC_N"] = N_p
    final["ATC_O"] = O_p
    final["ATC_S"] = S_p

    (final.round(2)).to_csv(out, index = None, encoding = 'utf-8')
########################################atom#################
def atc_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    atc('input_sam.csv',result_filename)
    os.remove('input_sam.csv')
def atc_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atc('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def atc_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atc('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def atc_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atc('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def atc_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atc('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def atc_st(file,out,N):
    readseq(file,'input_sam.csv')
    atom=pd.read_csv("atom.csv",header=None)
    #atom=pd.read_csv("atom.csv",header=None)
    at=pd.DataFrame()
    i = 0
    C_atom = []
    H_atom = []
    N_atom = []
    O_atom = []
    S_atom = []

    while i < len(atom):
        C_atom.append(atom.iloc[i,1].count("C"))
        H_atom.append(atom.iloc[i,1].count("H"))
        N_atom.append(atom.iloc[i,1].count("N"))
        O_atom.append(atom.iloc[i,1].count("O"))
        S_atom.append(atom.iloc[i,1].count("S"))
        i += 1
    atom["C_atom"]=C_atom
    atom["H_atom"]=H_atom
    atom["N_atom"]=N_atom
    atom["O_atom"]=O_atom
    atom["S_atom"]=S_atom

##############read file ##########
    data1 = pd.read_csv('input_sam.csv',header=None)
    data = pd.DataFrame(data1[0].str.upper())
    k1 = []
    for e in range(0,len(data)):
        s = 0
        k2 = []
        r = 0
        if len(data[0][e])%N == 0:
            k2.extend(repeat(int(len(data[0][e])/N),N))
        else:
            r = int(len(data[0][e])%N)
            k2.extend(repeat(int(len(data[0][e])/N),N-1))
            k2.append((int(len(data[0][e])/N))+r)
        for j in k2:
            df3 = data[0][e][s:j+s]
            k1.append(df3)
            s = j+s
    test = pd.DataFrame(k1)
    count_C = 0
    count_H = 0
    count_N = 0
    count_O = 0
    count_S = 0
    count = 0
    i1 = 0
    j = 0
    k = 0
    C_ct = []
    H_ct = []
    N_ct = []
    O_ct = []
    S_ct = []
    while i1 < len(test) :
        while j < len(test[0][i1]) :
            while k < len(atom) :
                if test.iloc[i1,0][j]==atom.iloc[k,0].replace(" ","") :
                    count_C = count_C + atom.iloc[k,2]
                    count_H = count_H + atom.iloc[k,3]
                    count_N = count_N + atom.iloc[k,4]
                    count_O = count_O + atom.iloc[k,5]
                    count_S = count_S + atom.iloc[k,6]
                k += 1
            k = 0
            j += 1
        C_ct.append(count_C)
        H_ct.append(count_H)
        N_ct.append(count_N)
        O_ct.append(count_O)
        S_ct.append(count_S)
        count_C = 0
        count_H = 0
        count_N = 0
        count_O = 0
        count_S = 0
        j = 0
        i1 += 1
    test["C_count"]=C_ct
    test["H_count"]=H_ct
    test["N_count"]=N_ct
    test["O_count"]=O_ct
    test["S_count"]=S_ct

    ct_total = []
    m = 0
    while m < len(test) :
        ct_total.append(test.iloc[m,1] + test.iloc[m,2] + test.iloc[m,3] + test.iloc[m,4] + test.iloc[m,5])
        m += 1
    test["count"]=ct_total
##########final output#####
    final = pd.DataFrame()
    n = 0
    p = 0
    C_p = []
    H_p = []
    N_p = []
    O_p = []
    S_p = []
    while n < len(test):
        C_p.append((test.iloc[n,1]/test.iloc[n,6])*100)
        H_p.append((test.iloc[n,2]/test.iloc[n,6])*100)
        N_p.append((test.iloc[n,3]/test.iloc[n,6])*100)
        O_p.append((test.iloc[n,4]/test.iloc[n,6])*100)
        S_p.append((test.iloc[n,5]/test.iloc[n,6])*100)
        n += 1
    final["C"] = C_p
    final["H"] = H_p
    final["N"] = N_p
    final["O"] = O_p
    final["S"] = S_p

    df3 = final
    bb = []
    for i in range(0,len(df3),N):
        aa = []
        for j in range(N):
            aa.extend(df3.loc[i+j])
        bb.append(aa)
    zz = pd.DataFrame(bb)
    header = []
    head = ['ATC_C_s','ATC_H_s','ATC_N_s','ATC_O_s','ATC_S_s']
    for e in range(1,N+1):
        for t in head:
            header.append(t+str(e))
    zz.columns=header
    (zz.round(2)).to_csv(out, index=None)
    os.remove('input_sam.csv')

#########################################bond#######################################################
def btc_wp(file,out):
    readseq(file,'input_sam.csv')
    tota = []
    hy = []
    Si = []
    Du = []
    b1 = []
    b2 = []
    b3 = []
    b4 = []
    bb = pd.DataFrame()
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv('input_sam.csv', header = None)
    bonds=pd.read_csv(paths_2+"/bonds.csv", sep = ",")
    for i in range(0,len(df)) :
        tot = 0
        h = 0
        S = 0
        D = 0
        tota.append([i])
        hy.append([i])
        Si.append([i])
        Du.append([i])
        for j in range(0,len(df[0][i])) :
            temp = df[0][i][j]
            for k in range(0,len(bonds)) :
                if bonds.iloc[:,0][k] == temp :
                    tot = tot + bonds.iloc[:,1][k]
                    h = h + bonds.iloc[:,2][k]
                    S = S + bonds.iloc[:,3][k]
                    D = D + bonds.iloc[:,4][k]
        tota[i].append(tot)
        hy[i].append(h)
        Si[i].append(S)
        Du[i].append(D)
    for m in range(0,len(df)) :
        b1.append(tota[m][1])
        b2.append(hy[m][1])
        b3.append(Si[m][1])
        b4.append(Du[m][1])

    bb["BTC_T"] = b1
    bb["BTC_H"] = b2
    bb["BTC_S"] = b3
    bb["BTC_D"] = b4

    bb.to_csv(out, index=None, encoding="utf-8")
    os.remove('input_sam.csv')
def btc_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    btc_wp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def btc_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    btc_wp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def btc_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    btc_wp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def btc_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    btc_wp('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def btc_st(file,out,n) :
    readseq(seq,'input_sam.csv')
    tota = []
    hy = []
    Si = []
    Du = []
    b1 = []
    b2 = []
    b3 = []
    b4 = []
    bb = pd.DataFrame()
    df = split('input_sam.csv',n)
    bonds=pd.read_csv(paths_2+"/bonds.csv")
    for i in range(0,len(df)) :
        tot = 0
        h = 0
        S = 0
        D = 0
        tota.append([i])
        hy.append([i])
        Si.append([i])
        Du.append([i])
        for j in range(0,len(df[0][i])) :
            temp = df[0][i][j]
            for k in range(0,len(bonds)) :
                if bonds.iloc[:,0][k] == temp :
                    tot = tot + bonds.iloc[:,1][k]
                    h = h + bonds.iloc[:,2][k]
                    S = S + bonds.iloc[:,3][k]
                    D = D + bonds.iloc[:,4][k]
        tota[i].append(tot)
        hy[i].append(h)
        Si[i].append(S)
        Du[i].append(D)
    for m in range(0,len(df)) :
        b1.append(tota[m][1])
        b2.append(hy[m][1])
        b3.append(Si[m][1])
        b4.append(Du[m][1])

    bb["total number of bonds"] = b1
    bb["hydrogen bonds"] = b2
    bb["single bonds"] = b3
    bb["double bonds"] = b4
    header = []
    header1 = ('BTC_T','BTC_H','BTC_S','BTC_D')
    for i in range(1,n+1):
        for j in header1:
            header.append(j+"_s"+str(i))
    qq = []
    for i in range(0,len(bb),n):
        aa = []
        for j in range(n):
            aa.extend(bb.loc[i+j])
        qq.append(aa)
    zz = pd.DataFrame(qq)
    zz.columns = header
    zz.to_csv(out, index=None)
    os.remove('input_sam.csv')
############################PhysicoChemical Properties###################################
PCP= pd.read_csv(paths_2+'/PhysicoChemical.csv', header=None)

headers = ['PCP_PC','PCP_NC','PCP_NE','PCP_PO','PCP_NP','PCP_AL','PCP_CY','PCP_AR','PCP_AC','PCP_BS','PCP_NE_pH','PCP_HB','PCP_HL','PCP_NT','PCP_HX','PCP_SC','PCP_SS_HE','PCP_SS_ST','PCP_SS_CO','PCP_SA_BU','PCP_SA_EX','PCP_SA_IN','PCP_TN','PCP_SM','PCP_LR','PCP_Z1','PCP_Z2','PCP_Z3','PCP_Z4','PCP_Z5'];

def encode(peptide):
    l=len(peptide);
    encoded=np.zeros(l);
    for i in range(l):
        if(peptide[i]=='A'):
            encoded[i] = 0;
        elif(peptide[i]=='C'):
            encoded[i] = 1;
        elif(peptide[i]=='D'):
            encoded[i] = 2;
        elif(peptide[i]=='E'):
            encoded[i] = 3;
        elif(peptide[i]=='F'):
            encoded[i] = 4;
        elif(peptide[i]=='G'):
            encoded[i] = 5;
        elif(peptide[i]=='H'):
            encoded[i] = 6;
        elif(peptide[i]=='I'):
            encoded[i] = 7;
        elif(peptide[i]=='K'):
            encoded[i] = 8;
        elif(peptide[i]=='L'):
            encoded[i] = 9;
        elif(peptide[i]=='M'):
            encoded[i] = 10;
        elif(peptide[i]=='N'):
            encoded[i] = 11;
        elif(peptide[i]=='P'):
            encoded[i] = 12;
        elif(peptide[i]=='Q'):
            encoded[i] = 13;
        elif(peptide[i]=='R'):
            encoded[i] = 14;
        elif(peptide[i]=='S'):
            encoded[i] = 15;
        elif(peptide[i]=='T'):
            encoded[i] = 16;
        elif(peptide[i]=='V'):
            encoded[i] = 17;
        elif(peptide[i]=='W'):
            encoded[i] = 18;
        elif(peptide[i]=='Y'):
            encoded[i] = 19;
        else:
            print('Wrong residue!');
    return encoded;
def lookup(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = encode(peptide);

    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);
def pcp_1(file,out123):

    if(type(file) == str):
        seq = pd.read_csv(file,header=None);
        #seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;

    l = len(seq);

    rows = PCP.shape[0]; # Number of features in our reference table
    col = 20 ; # Denotes the 20 amino acids

    seq=[seq[i].upper() for i in range(l)]
    sequenceFeature = [];
    sequenceFeature.append(headers); #To put property name in output csv

    for i in range(l): # Loop to iterate over each sequence
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): #Loop to iterate over each feature
            featureVal = lookup(seq[i],j)
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(round(featureVal/len(seq[i]),3));
            else:
                sequenceFeatureTemp.append('NaN')

        sequenceFeature.append(sequenceFeatureTemp);

    out = pd.DataFrame(sequenceFeature);
    file = open(out123,'w')
    with file:
        writer = csv.writer(file);
        writer.writerows(sequenceFeature);
    return sequenceFeature;

def pcp_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    pcp_1('input_sam.csv',result_filename)
    os.remove('input_sam.csv')
def pcp_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_1('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def pcp_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_1('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def pcp_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_1('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def pcp_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_1('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
###############################################################################################################################################
def pcp_st(file,out,v):
    readseq(file,'input_sam.csv')
    file1 = split('input_sam.csv',v)
    file2 = file1.iloc[:,0]
    pcp_1(file2,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*30)):
        pp.append(ss[i:i+(v*30)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
###############################RRI#################################
def RAAC(file,out):
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    count = 0
    cc = []
    i = 0
    x = 0
    temp = pd.DataFrame()
    f = open(out,'w')
    sys.stdout = f
    print("RRI_A,RRI_C,RRI_D,RRI_E,RRI_F,RRI_G,RRI_H,RRI_I,RRI_K,RRI_L,RRI_M,RRI_N,RRI_P,RRI_Q,RRI_R,RRI_S,RRI_T,RRI_V,RRI_W,RRI_Y,")
    for q in range(0,len(df1)):
        while i < len(std):
            cc = []
            for j in df1[0][q]:
                if j == std[i]:
                    count += 1
                    cc.append(count)
                else:
                    count = 0
            while x < len(cc) :
                if x+1 < len(cc) :
                    if cc[x]!=cc[x+1] :
                        if cc[x] < cc[x+1] :
                            cc[x]=0
                x += 1
            cc1 = [e for e in cc if e!= 0]
            cc = [e*e for e in cc if e != 0]
            zz= sum(cc)
            zz1 = sum(cc1)
            if zz1 != 0:
                zz2 = zz/zz1
            else:
                zz2 = 0
            print("%.2f"%zz2,end=',')
            i += 1
        i = 0
        print(" ")
    f.truncate()		

def rri_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    RAAC('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def rri_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    RAAC('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def rri_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    RAAC('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def rri_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    RAAC('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def rri_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    RAAC('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def rri_st(file,out,n):
    readseq(file,'input_sam.csv')
    filename, file_extension = os.path.splitext(file)
    file1 = split('input_sam.csv',n)
    data = file1
    count = 0
    cc = []
    i = 0
    x = 0
    temp = pd.DataFrame()
    f = open("sam.raac_split",'w')
    sys.stdout = f
    print("RRI_A,RRI_C,RRI_D,RRI_E,RRI_F,RRI_G,RRI_H,RRI_I,RRI_K,RRI_L,RRI_M,RRI_N,RRI_P,RRI_Q,RRI_R,RRI_S,RRI_T,RRI_V,RRI_W,RRI_Y,")
    for q in range(0,len(data)):
        mm = data[0][q]
        while i < len(std):
            cc = []
            for j in mm:
                if j == std[i]:
                    count += 1
                    cc.append(count)
                else:
                    count = 0
            while x < len(cc) :
                if x+1 < len(cc) :
                    if cc[x]!=cc[x+1] :
                        if cc[x] < cc[x+1] :
                            cc[x]=0
                x += 1
            cc1 = [e for e in cc if e!= 0]
            cc = [e*e for e in cc if e != 0]
            zz= sum(cc)
            zz1 = sum(cc1)
            if zz1 != 0:
                zz2 = zz/zz1
            else:
                zz2 = 0
            print("%.2f"%zz2,end=',')
            i += 1
        i = 0
        print(" ")
    f.truncate()
    df = pd.read_csv("sam.raac_split")
    df1 = df.iloc[:,:-1]
    header= []
    for i in range(1,n+1):
        for j in std:
            header.append('RRI_'+j+'_s'+str(i))
    ss = []
    for i in range(0,len(df1)):
        ss.extend(df1.loc[i])
    pp = []
    for i in range(0,len(ss),(n*20)):
        pp.append(ss[i:i+(n*20)])
    df5= pd.DataFrame(pp)
    df5.columns = header
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('sam.raac_split')
    os.remove('input_sam.csv')
##########################################PRI###########################
PCP= pd.read_csv(paths_2+'/PhysicoChemical.csv', header=None) #Our reference table for properties
headers_1 = ['PRI_PC','PRI_NC','PRI_NE','PRI_PO','PRI_NP','PRI_AL','PRI_CY','PRI_AR','PRI_AC','PRI_BS','PRI_NE_pH','PRI_HB','PRI_HL','PRI_NT','PRI_HX','PRI_SC','PRI_SS_HE','PRI_SS_ST','PRI_SS_CO','PRI_SA_BU','PRI_SA_EX','PRI_SA_IN','PRI_TN','PRI_SM','PRI_LR'];
def encode(peptide):
    l=len(peptide);
    encoded=np.zeros(l);
    for i in range(l):
        if(peptide[i]=='A'):
            encoded[i] = 0;
        elif(peptide[i]=='C'):
            encoded[i] = 1;
        elif(peptide[i]=='D'):
            encoded[i] = 2;
        elif(peptide[i]=='E'):
            encoded[i] = 3;
        elif(peptide[i]=='F'):
            encoded[i] = 4;
        elif(peptide[i]=='G'):
            encoded[i] = 5;
        elif(peptide[i]=='H'):
            encoded[i] = 6;
        elif(peptide[i]=='I'):
            encoded[i] = 7;
        elif(peptide[i]=='K'):
            encoded[i] = 8;
        elif(peptide[i]=='L'):
            encoded[i] = 9;
        elif(peptide[i]=='M'):
            encoded[i] = 10;
        elif(peptide[i]=='N'):
            encoded[i] = 11;
        elif(peptide[i]=='P'):
            encoded[i] = 12;
        elif(peptide[i]=='Q'):
            encoded[i] = 13;
        elif(peptide[i]=='R'):
            encoded[i] = 14;
        elif(peptide[i]=='S'):
            encoded[i] = 15;
        elif(peptide[i]=='T'):
            encoded[i] = 16;
        elif(peptide[i]=='V'):
            encoded[i] = 17;
        elif(peptide[i]=='W'):
            encoded[i] = 18;
        elif(peptide[i]=='Y'):
            encoded[i] = 19;
        else:
            print(peptide[i], ' is a wrong residue!');
    return encoded;
def lookup_1(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=[];
    peptide_num = encode(peptide);

    for i in range(l):
        out.append(PCP[peptide_num[i]][featureNum]);
    return out;
def binary_profile_1(file,featureNumb):
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;
    l = len(seq);
    bin_prof = [];
    for i in range(0,l):
        temp = lookup_1(seq[i],featureNumb);
        bin_prof.append(temp);
    return bin_prof;
def repeats(file,out123):
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        #seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;
    seq=[seq[i].upper() for i in range(len(seq))]
    dist =[];
    dist.append(headers_1);
    l = len(seq);
    for i in range(l):
        temp=[];
        for j in range(25):
            bin_prof = binary_profile_1(seq, j);
            if(j>=25):
                print('Error! Feature Number must be between 0-24');
                break;
            k=0;
            num=0;
            denom=0;
            ones=0;
            zeros=0;
            for j in range(len(bin_prof[i])):
                if(bin_prof[i][j]==0):
                    num+=k*k;
                    denom+=k;
                    k=0;
                    zeros+=1;
                elif(j==len(bin_prof[i])-1):
                    k+=1;
                    num+=k*k;
                    denom+=k;
                else:
                    k+=1;
                    ones+=1;
            if(ones!=0):
                answer = num/(ones*ones)
                temp.append(round(num/(ones*ones),2));
            elif(ones==0):
                temp.append(0);
        dist.append(temp)
    out = pd.DataFrame(dist)
    file1 = open(out123,'w')
    with file1:
        writer = csv.writer(file1);
        writer.writerows(dist);
    return out

def pri_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    repeats('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def pri_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    repeats('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def pri_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    repeats('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def pri_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    repeats('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def pri_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    repeats('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def rri_st(file,out,v):
    readseq(file,'input_sam.csv')
    file1 = split('input_sam.csv',v)
    file2 = file1.iloc[:,0]
    repeats(file2,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*25)):
        pp.append(ss[i:i+(v*25)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
	
#################################DDOR####################################
def DDOR(file,out) :
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    f = open(out,'w')
    sys.stdout = f
    for i in std:
        print('DDR_'+i, end=",")
    print("")
    for i in range(0,len(df1)):
        s = df1[0][i]
        p = s[::-1]
        for j in std:
            zz = ([pos for pos, char in enumerate(s) if char == j])
            pp = ([pos for pos, char in enumerate(p) if char == j])
            ss = []
            for i in range(0,(len(zz)-1)):
                ss.append(zz[i+1] - zz[i]-1)
            if zz == []:
                ss = []
            else:
                ss.insert(0,zz[0])
                ss.insert(len(ss),pp[0])
            cc1=  (sum([e for e in ss])+1)
            cc = sum([e*e for e in ss])
            zz2 = cc/cc1
            print("%.2f"%zz2,end=",")
        print("")
    f.truncate()
def ddr_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    DDOR('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ddr_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    DDOR('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ddr_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    DDOR('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ddr_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    DDOR('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ddr_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    DDOR('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.iloc[:,:-1].to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ddr_st(file,out,v):
    readseq(file,'input_sam.csv')
    file1 = split('input_sam.csv',v)
    df = file1
    df1 = pd.DataFrame(df[0].str.upper())
    f = open('tempfile_out','w')
    sys.stdout = f
    for i in std:
        print('DDR_'+i, end=",")
    print("")
    for i in range(0,len(df1)):
        s = df1[0][i]
        p = s[::-1]
        for j in std:
            zz = ([pos for pos, char in enumerate(s) if char == j])
            pp = ([pos for pos, char in enumerate(p) if char == j])
            ss = []
            for i in range(0,(len(zz)-1)):
                ss.append(zz[i+1] - zz[i]-1)
            if zz == []:
                ss = []
            else:
                ss.insert(0,zz[0])
                ss.insert(len(ss),pp[0])
            cc1=  (sum([e for e in ss])+1)
            cc = sum([e*e for e in ss])
            zz2 = cc/cc1
            print("%.2f"%zz2,end=",")
        print("")
    f.truncate()
    df3 = pd.read_csv("tempfile_out")
    df4 = df3.iloc[:,:-1]
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*20)):
        pp.append(ss[i:i+(v*20)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')

########################Shannon_Entropy Whole protein######################################
def entropy_single(seq):
    seq=seq.upper()
    num, length = Counter(seq), len(seq)
    return -sum( freq/length * math.log(freq/length, 2) for freq in num.values())

def SE(filename,out):
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    Val=[]
    header=["SEP"]
    for i in range(len(data)):
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Error: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        Val.append(round((entropy_single(str(data[i]))),3))
        #print(Val[i])
        file= open(out,'w', newline='\n')#output file
        with file:
            writer=csv.writer(file,delimiter='\n');
            writer.writerow(header)
            writer.writerow(Val);
    return Val

def sep_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    SE('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def sep_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def sep_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def sep_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def sep_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def sep_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE_split('sam_input.csv',sp,result_filename)
    os.remove('input_sam.csv')

def SE_split(file,v,out):
    SE(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),v):
        pp.append(ss[i:i+v])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

################################Shannon_Entropy residue#############################################
def SE_residue_level(filename,out):
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    data2=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    Val=np.zeros(len(data))
    GH=[]
    for i in range(len(data)):
        my_list={'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Error: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        seq=data[i]
        seq=seq.upper()
        num, length = Counter(seq), len(seq)
        num=dict(sorted(num.items()))
        C=list(num.keys())
        F=list(num.values())
        for key, value in my_list.items():
             for j in range(len(C)):
                if key == C[j]:
                    my_list[key] = round(((F[j]/length)* math.log(F[j]/length, 2)),3)
        GH.append(list(my_list.values()))
    file= open(out,'w', newline='')#output file
    with file:
        writer=csv.writer(file);
        writer.writerow(('SER_A','SER_C','SER_D','SER_E','SER_F','SER_G','SER_H','SER_I','SER_K','SER_L','SER_M','SER_N','SER_P','SER_Q','SER_R','SER_S','SER_T','SER_V','SER_W','SER_Y'));
        writer.writerows(GH);
    return(GH)

def SE_residue_level_split(file,v,out):
    SE_residue_level(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*20)):
        pp.append(ss[i:i+(v*20)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

def ser_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    SE_residue_level('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ser_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE_residue_level('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ser_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE_residue_level('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ser_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE_residue_level('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ser_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE_residue_level('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ser_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    SE_residue_level_split('sam_input.csv',sp,result_filename)
    os.remove('input_sam.csv')

#########################Shanon entropy for PCP#################################
def lookup(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = encode(peptide);
    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);
def pcp(file):
    SEP_headers = ['SEP_PC','SEP_NC','SEP_NE','SEP_PO','SEP_NP','SEP_AL','SEP_CY','SEP_AR','SEP_AC','SEP_BS','SEP_NE_pH','SEP_HB','SEP_HL','SEP_NT','SEP_HX','SEP_SC','SEP_SS_HE','SEP_SS_ST','SEP_SS_CO','SEP_SA_BU','SEP_SA_EX','SEP_SA_IN','SEP_TN','SEP_SM','SEP_LR']
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;
    l = len(seq);
    rows = PCP.shape[0]; # Number of features in our reference table
    col = 20 ; # Denotes the 20 amino acids
    seq=[seq[i].upper() for i in range(l)]
    sequenceFeature = [];
    sequenceFeature.append(SEP_headers); #To put property name in output csv
    
    for i in range(l): # Loop to iterate over each sequence
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): #Loop to iterate over each feature
            featureVal = lookup(seq[i],j)   
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(featureVal/len(seq[i]))
            else:
                sequenceFeatureTemp.append('NaN')
        sequenceFeature.append(sequenceFeatureTemp);
    out = pd.DataFrame(sequenceFeature);
    return sequenceFeature;
def phyChem(file,mode='all',m=0,n=0):
    if(type(file) == str):
        seq1 = pd.read_csv(file,header=None, sep=',');
        seq1 = pd.DataFrame(seq1[0].str.upper())
        seq=[]
        [seq.append(seq1.iloc[i][0]) for i in range(len(seq1))]
    else:
        seq  = file;
    l = len(seq);
    newseq = [""]*l; # To store the n-terminal sequence
    for i in range(0,l):
        l = len(seq[i]);
        if(mode=='NT'):
            n=m;
            if(n!=0):
                newseq[i] = seq[i][0:n];
            elif(n>l):
                print('Warning! Sequence',i,"'s size is less than n. The output table would have NaN for this sequence");
            else:
                print('Value of n is mandatory, it cannot be 0')
                break;
        elif(mode=='CT'):
            n=m;
            if(n!=0):
                newseq[i] = seq[i][(len(seq[i])-n):]
            elif(n>l):
                print('WARNING: Sequence',i+1,"'s size is less than the value of n given. The output table would have NaN for this sequence");
            else:
                print('Value of n is mandatory, it cannot be 0')
                break;
        elif(mode=='all'):
            newseq = seq;
        elif(mode=='rest'):
            if(m==0):
                print('Kindly provide start index for rest, it cannot be 0');
                break;
            else:
                if(n<=len(seq[i])):
                    newseq[i] = seq[i][m-1:n+1]
                elif(n>len(seq[i])):
                    newseq[i] = seq[i][m-1:len(seq[i])]
                    print('WARNING: Since input value of n for sequence',i+1,'is greater than length of the protein, entire sequence starting from m has been considered')
        else:
            print("Wrong Mode. Enter 'NT', 'CT','all' or 'rest'");        
    output = pcp(newseq);
    return output
def shannons(filename,out123):
    SEP_headers = ['SEP_PC','SEP_NC','SEP_NE','SEP_PO','SEP_NP','SEP_AL','SEP_CY','SEP_AR','SEP_AC','SEP_BS','SEP_NE_pH','SEP_HB','SEP_HL','SEP_NT','SEP_HX','SEP_SC','SEP_SS_HE','SEP_SS_ST','SEP_SS_CO','SEP_SA_BU','SEP_SA_EX','SEP_SA_IN','SEP_TN','SEP_SM','SEP_LR']
    if(type(filename) == str):
        seq1 = pd.read_csv(filename,header=None, sep=',');
        seq1 = pd.DataFrame(seq1[0].str.upper())
    else:
        seq1  = filename;
    seq=[]
    [seq.append(seq1.iloc[i][0]) for i in range(len(seq1))]
    comp = phyChem(seq);
    new = [comp[i][0:25] for i in range(len(comp))]
    entropy  = [];
    entropy.append(SEP_headers[0:25])
    for i in range(1,len(new)):
        seqEntropy = [];
        for j in range(len(new[i])):
            p = new[i][j]; 
            if((1-p) == 0. or p==0.):
                temp = 0;#to store entropy of each sequence
            else:
                temp = -(p*math.log2(p)+(1-p)*math.log2(1-p));
            seqEntropy.append(round(temp,3));
        entropy.append(seqEntropy);
    out = pd.DataFrame(entropy);
    file = open(out123,'w')
    with file:
        writer = csv.writer(file);
        writer.writerows(entropy);
    return entropy;
def shannons_split(file,v,out):
    shannons(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*25)):
        pp.append(ss[i:i+(v*25)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')	
	
def spc_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    shannons('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df = df.round(3)
    df.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def spc_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    shannons('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def spc_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    shannons('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def spc_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    shannons('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def spc_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    shannons('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def spc_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    shannons_split('sam_input.csv',sp,result_filename)
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
####################autocorr####################
def p_aa(prop,a):
    if ((a.upper()=='A') or (a.upper()=='C') or (a.upper()=='D') or (a.upper()=='E') or (a.upper()=='F') or (a.upper()=='G') or (a.upper()=='H') or (a.upper()=='I') or (a.upper()=='K') or (a.upper()=='L') or (a.upper()=='M') or (a.upper()=='N') or (a.upper()=='P') or (a.upper()=='Q') or (a.upper()=='R') or (a.upper()=='S') or (a.upper()=='T') or (a.upper()=='V') or (a.upper()=='W') or (a.upper()=='Y')):
        data=pd.read_table(paths_2+'/z_aaindex.csv',sep=',',index_col='INDEX' )
        p=data.loc[prop][a.upper()]
        return p
    else:
        print("Error!: Invalid sequence. Special character(s)/invalid amino acid letter(s) present in input.")
        return
def NMB(prop,seq,d):
    if (d<=30):
        sum=0
        for i in range(len(seq)-d):
            sum=sum+p_aa(prop,seq[i])*p_aa(prop,seq[i+d])
        ats=sum/(len(seq)-d)
        return ats
    else:
        print("Error!: d should be less than equal to 30")
        return
def pav(prop,seq):
    av=0
    for i in range(len(seq)):
        av=av+p_aa(prop,seq[i])
    av=av/len(seq)
    return av
def moran(prop,seq,d):
    if (d<=30):
        s1=0
        s2=0
        p_bar=pav(prop,seq)
        for i in range(len(seq)-d):
            s1=s1+(p_aa(prop,seq[i])-p_bar)*(p_aa(prop,seq[i+d])-p_bar)
        for i in range(len(seq)):
            s2=s2+(p_aa(prop,seq[i])-p_bar)*(p_aa(prop,seq[i])-p_bar)
        return (s1/(len(seq)-d))/(s2/len(seq))
    else:
        print("Error!: d should be less than equal to 30")
        return
def geary(prop,seq,d):
    if (d<=30):
        s1=0
        s2=0
        p_bar=pav(prop,seq)
        for i in range(len(seq)-d):
            s1=s1+(p_aa(prop,seq[i])-p_aa(prop,seq[i+d]))*(p_aa(prop,seq[i])-p_aa(prop,seq[i+d]))
        for i in range(len(seq)):
            s2=s2+(p_aa(prop,seq[i])-p_bar)*(p_aa(prop,seq[i])-p_bar)
        return (s1/(2*(len(seq)-d)))/(s2/(len(seq)-1))
    else:
        print("Error!: d should be less than equal to 30")
        return

def autocorr_full_aa(filename,d,out):
    if (d<=30):
        seq_data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
        prop=list((pd.read_csv(paths_2+'/aaindices.csv',sep=',',header=None)).iloc[0,:])
        output=[[]]
        for k in range(len(prop)):
            output[0]=output[0]+['ACR'+str(d)+'_'+'MB','ACR'+str(d)+'_'+'MO','ACR'+str(d)+'_'+'GE']
        temp=[]
        for i in range(len(seq_data)):
            for j in range(len(prop)):
                temp=temp+[round(NMB(prop[j],seq_data[i],d),3),round(moran(prop[j],seq_data[i],d),3),round(geary(prop[j],seq_data[i],d),3)]
            output.append(temp)
            temp=[]
        file = open(out,'w')
        with file:
            writer = csv.writer(file);
            writer.writerows(output);
        return output
    else:
        print("Error!: d should be less than equal to 30")
        return

def autocorr_split(file,v,lg,out):
    autocorr_full_aa(file,lg,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*3)):
        pp.append(ss[i:i+(v*3)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

def acr_wp(seq,result_filename,lg):
    readseq(seq,'input_sam.csv')
    autocorr_full_aa('input_sam.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def acr_nt(seq,result_filename,n,lg):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def acr_ct(seq,result_filename,c,lg):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def acr_rt(seq,result_filename,n,c,lg):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def acr_nct(seq,result_filename,n,lg):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    autocorr_full_aa('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def acr_st(seq,result_filename,sp,lg):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    autocorr_split('sam_input.csv',sp,lg,result_filename)
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')

##################paac####################
def val(AA_1, AA_2, aa, mat):
    return sum([(mat[i][aa[AA_1]] - mat[i][aa[AA_2]]) ** 2 for i in range(len(mat))]) / len(mat)
def paac_1(file,lambdaval,w=0.05):
    data1 = pd.read_csv(paths_2+"/data", sep = "\t")
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    dd = []
    cc = []
    pseudo = []
    aa = {}
    for i in range(len(std)):
        aa[std[i]] = i
    for i in range(0,3):
        mean = sum(data1.iloc[i][1:])/20
        rr = math.sqrt(sum([(p-mean)**2 for p in data1.iloc[i][1:]])/20)
        dd.append([(p-mean)/rr for p in data1.iloc[i][1:]])
        zz = pd.DataFrame(dd)
    head = []
    for n in range(1, lambdaval + 1):
        head.append('_lam' + str(n))
    head = ['PAAC'+str(lambdaval)+sam for sam in head]
    pp = pd.DataFrame()
    ee = []
    for k in range(0,len(df1)):
        cc = []
        pseudo1 = [] 
        for n in range(1,lambdaval+1):
            cc.append(sum([val(df1[0][k][p], df1[0][k][p + n], aa, dd) for p in range(len(df1[0][k]) - n)]) / (len(df1[0][k]) - n))
            qq = pd.DataFrame(cc)
        pseudo = pseudo1 + [(w * p) / (1 + w * sum(cc)) for p in cc]
        ee.append(pseudo)
        ii = round(pd.DataFrame(ee, columns = head),4)
    ii.to_csv(filename+".lam",index = None)
		
def paac(file,lambdaval,out,w=0.05):
    filename, file_extension = os.path.splitext(file)
    paac_1(file,lambdaval,w=0.05)
    aac_comp(file,filename+".aac")
    data1 = pd.read_csv(filename+".aac")
    header = ['PAAC'+str(lambdaval)+'_A','PAAC'+str(lambdaval)+'_C','PAAC'+str(lambdaval)+'_D','PAAC'+str(lambdaval)+'_E','PAAC'+str(lambdaval)+'_F','PAAC'+str(lambdaval)+'_G','PAAC'+str(lambdaval)+'_H','PAAC'+str(lambdaval)+'_I','PAAC'+str(lambdaval)+'_K','PAAC'+str(lambdaval)+'_L','PAAC'+str(lambdaval)+'_M','PAAC'+str(lambdaval)+'_N','PAAC'+str(lambdaval)+'_P','PAAC'+str(lambdaval)+'_Q','PAAC'+str(lambdaval)+'_R','PAAC'+str(lambdaval)+'_S','PAAC'+str(lambdaval)+'_T','PAAC'+str(lambdaval)+'_V','PAAC'+str(lambdaval)+'_W','PAAC'+str(lambdaval)+'_Y','Un']	
    data1.columns = header    
    data2 = pd.read_csv(filename+".lam")
    data3 = pd.concat([data1.iloc[:,:-1],data2], axis = 1).reset_index(drop=True)
    data3.to_csv(out,index=None)
    os.remove(filename+".lam")
    os.remove(filename+".aac")


def paac_wp(seq,result_filename,lg,pw):
    readseq(seq,'input_sam.csv')
    paac('input_sam.csv',lg,result_filename,pw)
    os.remove('input_sam.csv')
def paac_nt(seq,result_filename,n,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    paac('sam_input.csv',lg,'tempfile_out',pw)
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def paac_ct(seq,result_filename,c,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    paac('sam_input.csv',lg,'tempfile_out',pw)
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def paac_rt(seq,result_filename,n,c,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    paac('sam_input.csv',lg,'tempfile_out',pw)
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def paac_nct(seq,result_filename,n,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    paac('sam_input.csv',lg,'tempfile_out',pw)
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def paac_st(seq,result_filename,sp,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    paac_split('sam_input.csv',sp,lg,result_filename,pw)
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
def paac_split(file,v,lg,out,w):
    paac(file,lg,'tempfile_out',w)
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')
	
######################apaac############################	
def apaac_1(file,lambdaval,w=0.05):
    data1 = pd.read_csv(paths_2+"/data", sep = "\t")
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    dd = []
    cc = []
    pseudo = []
    aa = {}
    for i in range(len(std)):
        aa[std[i]] = i
    for i in range(0,3):
        mean = sum(data1.iloc[i][1:])/20
        rr = math.sqrt(sum([(p-mean)**2 for p in data1.iloc[i][1:]])/20)
        dd.append([(p-mean)/rr for p in data1.iloc[i][1:]])
        zz = pd.DataFrame(dd)
    head = []
    for n in range(1, lambdaval + 1):
        for e in ('HB','HL','SC'):
            head.append(e+'_lam' + str(n))
    head = ['APAAC'+str(lambdaval)+'_'+sam for sam in head]
    pp = pd.DataFrame()
    ee = []
    for k in range(0,len(df1)):
        cc = [] 
        for n in range(1,lambdaval+1):
            for b in range(0,len(zz)):
                cc.append(sum([zz.loc[b][aa[df1[0][k][p]]] * zz.loc[b][aa[df1[0][k][p + n]]] for p in range(len(df1[0][k]) - n)]) / (len(df1[0][k]) - n))
                qq = pd.DataFrame(cc)
        pseudo = [(w * p) / (1 + w * sum(cc)) for p in cc]
        ee.append(pseudo)
        ii = round(pd.DataFrame(ee, columns = head),4)
    ii.to_csv(filename+".plam",index = None)
		
def apaac(file,lambdaval,out,w=0.05):
    filename, file_extension = os.path.splitext(file)
    apaac_1(file,lambdaval,w=0.05)
    aac_comp(file,filename+".aac")
    data1 = pd.read_csv(filename+".aac")
    headaac = []
    for i in std:
        headaac.append('APAAC'+str(lambdaval)+'_'+i)
    headaac.insert(len(headaac),0)
    data1.columns = headaac
    data2 = pd.read_csv(filename+".plam")
    data3 = pd.concat([data1.iloc[:,:-1],data2], axis = 1).reset_index(drop=True)
    data3.to_csv(out, index = None)
    os.remove(filename+".plam")
    os.remove(filename+".aac")
def apaac_split(file,v,lg,out,w):
    apaac(file,lg,'tempfile_out',w)
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

def apaac_wp(seq,result_filename,lg,pw):
    readseq(seq,'input_sam.csv')
    apaac('input_sam.csv',lg,result_filename,pw)
    os.remove('input_sam.csv')
def apaac_nt(seq,result_filename,n,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    apaac('sam_input.csv',lg,'tempfile_out',pw)
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def apaac_ct(seq,result_filename,c,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    apaac('sam_input.csv',lg,'tempfile_out',pw)
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def apaac_rt(seq,result_filename,n,c,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    apaac('sam_input.csv',lg,'tempfile_out',pw)
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def apaac_nct(seq,result_filename,n,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    apaac('sam_input.csv',lg,'tempfile_out',pw)
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def apaac_st(seq,result_filename,sp,lg,pw):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    apaac_split('sam_input.csv',sp,lg,result_filename,pw)
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
###################################qos#######################################
def qos(file,gap,out,w=0.1):
    ff = []
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df[0].str.upper())
    for i in range(0,len(df2)):
        ff.append(len(df2[0][i]))
    if min(ff) < gap:
        print("Error: All sequences' length should be higher than :", gap)
    else:
        mat1 = pd.read_csv(paths_2+"/Schneider-Wrede.csv", index_col = 'Name')
        mat2 = pd.read_csv(paths_2+"/Grantham.csv", index_col = 'Name')
        s1 = []
        s2 = []
        for i in range(0,len(df2)):
            for n in range(1, gap+1):
                sum1 = 0
                sum2 = 0
                for j in range(0,(len(df2[0][i])-n)):
                    sum1 = sum1 + (mat1[df2[0][i][j]][df2[0][i][j+n]])**2
                    sum2 = sum2 + (mat2[df2[0][i][j]][df2[0][i][j+n]])**2
                s1.append(sum1)
                s2.append(sum2)
        zz = pd.DataFrame(np.array(s1).reshape(len(df2),gap))
        zz["sum"] = zz.sum(axis=1)
        zz2 = pd.DataFrame(np.array(s2).reshape(len(df2),gap))
        zz2["sum"] = zz2.sum(axis=1)
        c1 = []
        c2 = []
        c3 = []
        c4 = []
        h1 = []
        h2 = []
        h3 = []
        h4 = []
        for aa in std:
            h1.append('QSO'+str(gap)+'_SC_' + aa)
        for aa in std:
            h2.append('QSO'+str(gap)+'_G_' + aa)
        for n in range(1, gap+1):
            h3.append('SC' + str(n))
        h3 = ['QSO'+str(gap)+'_'+sam for sam in h3]
        for n in range(1, gap+1):
            h4.append('G' + str(n))
        h4 = ['QSO'+str(gap)+'_'+sam for sam in h4]
        for i in range(0,len(df2)):
            AA = {}
            for j in std:
                AA[j] = df2[0][i].count(j)
                c1.append(AA[j] / (1 + w * zz['sum'][i]))
                c2.append(AA[j] / (1 + w * zz2['sum'][i]))
            for k in range(0,gap):
                c3.append((w * zz[k][i]) / (1 + w * zz['sum'][i]))
                c4.append((w * zz[k][i]) / (1 + w * zz['sum'][i]))
        pp1 = np.array(c1).reshape(len(df2),len(std))
        pp2 = np.array(c2).reshape(len(df2),len(std))
        pp3 = np.array(c3).reshape(len(df2),gap)
        pp4 = np.array(c4).reshape(len(df2),gap)
        zz5 = round(pd.concat([pd.DataFrame(pp1, columns = h1),pd.DataFrame(pp2,columns = h2),pd.DataFrame(pp3, columns = h3),pd.DataFrame(pp4, columns = h4)], axis = 1),4)
        zz5.to_csv(out, index = None, encoding = 'utf-8')	
def qos_split(file,v,lg,out,w):
    qos(file,lg,'tempfile_out',w)
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

def qos_wp(seq,result_filename,lg,wq):
    readseq(seq,'input_sam.csv')
    qos('input_sam.csv',lg,result_filename,wq)
    os.remove('input_sam.csv')
def qos_nt(seq,result_filename,n,lg,wq):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    qos('sam_input.csv',lg,'tempfile_out',wq)
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def qos_ct(seq,result_filename,c,lg,wq):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    qos('sam_input.csv',lg,'tempfile_out',wq)
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def qos_rt(seq,result_filename,n,c,lg,wq):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    qos('sam_input.csv',lg,'tempfile_out',wq)
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def qos_nct(seq,result_filename,n,lg,wq):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    qos('sam_input.csv',lg,'tempfile_out',wq)
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def qos_st(seq,result_filename,sp,lg,wq):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    qos_split('sam_input.csv',sp,lg,result_filename,wq)
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
##########################soc################
def soc(file,gap,out):
    ff = []
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df[0].str.upper())
    for i in range(0,len(df2)):
        ff.append(len(df2[0][i]))
    if min(ff) < gap:
        print("Error: All sequences' length should be higher than :", gap)
        return 0
    mat1 = pd.read_csv(paths_2+"/Schneider-Wrede.csv", index_col = 'Name')
    mat2 = pd.read_csv(paths_2+"/Grantham.csv", index_col = 'Name')
    h1 = []
    h2 = []
    for n in range(1, gap+1):
        h1.append('SC' + str(n))
    for n in range(1, gap + 1):
        h2.append('G' + str(n))
    h1 = ['SOC'+str(gap)+'_'+sam for sam in h1]
    h2 = ['SOC'+str(gap)+'_'+sam for sam in h2]
    s1 = []
    s2 = []
    for i in range(0,len(df2)):
        for n in range(1, gap+1):
            sum = 0
            sum1 =0
            sum2 =0
            sum3 =0
            for j in range(0,(len(df2[0][i])-n)):
                sum = sum + (mat1[df2[0][i][j]][df2[0][i][j+n]])**2
                sum1 = sum/(len(df2[0][i])-n)
                sum2 = sum2 + (mat2[df2[0][i][j]][df2[0][i][j+n]])**2
                sum3 = sum2/(len(df2[0][i])-n)
            s1.append(sum1)
            s2.append(sum3)
    zz = np.array(s1).reshape(len(df2),gap)
    zz2 = np.array(s2).reshape(len(df2),gap)
    zz3 = round(pd.concat([pd.DataFrame(zz, columns = h1),pd.DataFrame(zz2,columns = h2)], axis = 1),4)
    zz3.to_csv(out, index = None, encoding = 'utf-8') 

def soc_split(file,v,lg,out):
    soc(file,lg,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

def soc_wp(seq,result_filename,lg):
    readseq(seq,'input_sam.csv')
    soc('input_sam.csv',lg,result_filename)
    os.remove('input_sam.csv')
def soc_nt(seq,result_filename,n,lg):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    soc('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def soc_ct(seq,result_filename,c,lg):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    soc('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def soc_rt(seq,result_filename,n,c,lg):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    soc('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def soc_nct(seq,result_filename,n,lg):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    soc('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def soc_st(seq,result_filename,sp,lg):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    soc_split('sam_input.csv',sp,lg,result_filename)
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
##########################################CTC###################################
x = [1, 2, 3, 4, 5, 6,7]
p=[]
Y=[]
LS=[]


for i in range(len(x)):
    p=itertools.product(x,repeat=3)
    p=list(p)

def concatenate_list_data(list):
    result= ''
    for element in list:
        result += str(element)
    return result

for i in range(len(p)):
    LS.append(concatenate_list_data(p[i]))

def repstring(string):
    string=string.upper()
    char={"A":"1","G":"1","V":"1","I":"2","L":"2","F":"2","P":"2","Y":"3","M":"3","T":"3","S":"3","H":"4","N":"4","Q":"4","W":"4","R":"5","K":"5","D":"6","E":"6","C":"7"}
    string=list(string)
    for index,item in enumerate(string):
        for key,value in char.items():
            if item==key:
                string[index]=value
    return("".join(string))

def occurrences(string, sub_string):
    count=0
    beg=0
    while(string.find(sub_string,beg)!=-1) :
        count=count+1
        beg=string.find(sub_string,beg)
        beg=beg+1
    return count


def CTC(filename,out):
    df = pd.DataFrame(columns=['Sequence','Triad:Frequency'])
    data=list((pd.read_csv(filename,sep=',',header=None)).iloc[:,0])
    for i in range(len(data)):
        data1=''
        data1=str(data[i])
        data1=data1.upper()
        allowed = set(('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
        is_data_invalid = set(data1).issubset(allowed)
        if is_data_invalid==False:
            print("Errror: Please check for invalid inputs in the sequence.","\nError in: ","Sequence number=",i+1,",","Sequence = ",data[i],",","\nNOTE: Spaces, Special characters('[@_!#$%^&*()<>?/\|}{~:]') and Extra characters(BJOUXZ) should not be there.")
            return
        df.at[i,'Sequence'] = data[i]
        Y.append("".join(repstring(str(data[i]))))
    val2=[[]]
    for f in range(len(LS)):
        val2[0]=val2[0]+["CTC_"+str(LS[f])]
    for j in range(len(data)):
        MM=[]
        for m in range(len(LS)):
            MM=MM+[occurrences(Y[j],LS[m])]
        Min_MM=min(MM)
        Max_MM=max(MM)
        if (Max_MM==0):
            print("Errror: Splits/ Sequence length should be greater than equal to 3")
            return
        val=[]
#         val.append(data[j])
        for k in range(len(LS)):
            val=val+[round(((occurrences(Y[j],LS[k])-Min_MM)/Max_MM),3)]
        val2.append(val)
#     print(val2)
    #file= open(sys.argv[2],'w', newline='')#output file
    file= open(out,'w', newline='')
    with file:
        writer=csv.writer(file);
        writer.writerows(val2);
    return val2

def ctc_split(file,v,out):
    CTC(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

def ctc_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    CTC('input_sam.csv',result_filename)
    os.remove('input_sam.csv')
def ctc_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    CTC('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ctc_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    CTC('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ctc_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    CTC('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ctc_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    CTC('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ctc_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    ctc_split('sam_input.csv',sp,result_filename)
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
######################################AAC###################################
def ctd(file,out):
    attr=pd.read_csv(paths_2+"/aa_attr_group.csv", sep="\t")
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df = pd.DataFrame(df1[0].str.upper())
    n = 0
    stt1 = []
    m = 1
    for i in range(0,len(attr)) :
        st =[]
        stt1.append([])
        for j in range(0,len(df)) :
            stt1[i].append([])
            for k in range(0,len(df[0][j])) :
                while m < 4 :
                    while n < len(attr.iloc[i,m]) :
                        if df[0][j][k] == attr.iloc[i,m][n] :
                            st.append(m)
                            stt1[i][j].append(m)
                        n += 2
                    n = 0
                    m += 1
                m = 1
#####################Composition######################
    f = open("compout_1", 'w')
    sys.stdout = f
    std = [1,2,3]
    print("1,2,3,")
    for p in range (0,len(df)) :
        for ii in range(0,len(stt1)) :
            #for jj in stt1[ii][p]:
            for pp in std :
                count = 0
                for kk in stt1[ii][p] :
                    temp1 = kk
                    if temp1 == pp :
                        count += 1
                    composition = (count/len(stt1[ii][p]))*100
                print("%.2f"%composition, end = ",")
            print("")
    f.truncate()

#################################Transition#############
    tt = []
    tr=[]
    kk =0
    for ii in range(0,len(stt1)) :
        tt = []
        tr.append([])
        for p in range (0,len(df)) :
            tr[ii].append([])
            while kk < len(stt1[ii][p]) :
                if kk+1 <len(stt1[ii][p]):
                #if  stt1[ii][p][kk] < stt1[ii][p][kk+1] or stt1[ii][p][kk] > stt1[ii][p][kk+1]: # condition for adjacent values
                    tt.append(stt1[ii][p][kk])
                    tt.append(stt1[ii][p][kk+1])
                    tr[ii][p].append(stt1[ii][p][kk])
                    tr[ii][p].append(stt1[ii][p][kk+1])

                kk += 1
            kk = 0

    pp = 0
    xx = []
    xxx = []
    for mm in range(0,len(tr)) :
        xx = []
        xxx.append([])
        for nn in range(0,len(tr[mm])):
            xxx[mm].append([])
            while pp < len(tr[mm][nn]) :
                xx .append(tr[mm][nn][pp:pp+2])
                xxx[mm][nn].append(tr[mm][nn][pp:pp+2])
                pp+=2
            pp = 0

    f1 = open("compout_2", 'w')
    sys.stdout = f1
    std1 = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
    print("1->1,1->2,1->3,2->1,2->2,2->3,3->1,3->2,3->3,")
    for rr in range(0,len(df)) :
        for qq in range(0,len(xxx)):
            for tt in std1 :
                count = 0
                for ss in xxx[qq][rr] :
                    temp2 = ss
                    if temp2 == tt :
                        count += 1
                print(count, end = ",")
            print("")
    f1.truncate()

    #################################Distribution#############
    c_11 = []
    c_22 = []
    c_33 = []
    zz = []
    #print("0% 25% 50% 75% 100%")
    for x in range(0,len(stt1)) :
        #c_11.append([])
        c_22.append([])
        #c_33.append([])
        yy_c_1 = []
        yy_c_2 = []
        yy_c_3 = []
        ccc = []

        k = 0
        j = 0
        for y in range(0,len(stt1[x])):
            #c_11[x].append([])
            c_22[x].append([])
            for i in range(1,4) :
                cc = []
                c1 = [index for index,value in enumerate(stt1[x][y]) if value == i]
                c_22[x][y].append(c1)
    cc = []
    for ss in range(0,len(df)):
        for uu in range(0,len(c_22)):
            for mm in range(0,3):
                for ee in range(0,101,25):
                    k = (ee*(len(c_22[uu][ss][mm])))/100
                    cc.append(math.floor(k))
    f2 = open('compout_3', 'w')
    sys.stdout = f2
    print("0% 25% 50% 75% 100%")
    for i in range (0,len(cc),5):
        print(*cc[i:i+5])
    f2.truncate()
    head = []
    header1 = ['CeTD_HB','CeTD_VW','CeTD_PO','CeTD_PZ','CeTD_CH','CeTD_SS','CeTD_SA']
    for i in header1:
        for j in range(1,4):
            head.append(i+str(j))
    df11 = pd.read_csv("compout_1")
    df_1 = df11.iloc[:,:-1]
    zz = pd.DataFrame()
    for i in range(0,len(df_1),7):
        zz = zz.append(pd.DataFrame(pd.concat([df_1.loc[i],df_1.loc[i+1],df_1.loc[i+2],df_1.loc[i+3],df_1.loc[i+4],df_1.loc[i+5],df_1.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
    zz.columns = head
    #zz.to_csv(filename+".ctd_comp", index=None, encoding='utf-8')
    head2 = []
    header2 = ['CeTD_11','CeTD_12','CeTD_1-3','CeTD_21','CeTD_22','CeTD_23','CeTD_31','CeTD_32','CeTD_33']
    for i in header2:
        for j in ('HB','VW','PO','PZ','CH','SS','SA'):
            head2.append(i+'_'+str(j))
    df12 = pd.read_csv("compout_2")
    df_2 = df12.iloc[:,:-1]
    ss = pd.DataFrame()
    for i in range(0,len(df_2),7):
        ss = ss.append(pd.DataFrame(pd.concat([df_2.loc[i],df_2.loc[i+1],df_2.loc[i+2],df_2.loc[i+3],df_2.loc[i+4],df_2.loc[i+5],df_2.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
    ss.columns = head2
    #ss.to_csv(filename+".ctd_trans", index=None, encoding='utf-8')
    head3 = []
    header3 = ['CeTD_0_p','CeTD_25_p','CeTD_50_p','CeTD_75_p','CeTD_100_p']
    header4 = ['HB','VW','PO','PZ','CH','SS','SA']
    for j in range(1,4):
        for k in header4:
            for i in header3:
                head3.append(i+'_'+k+str(j))
    df_3 = pd.read_csv("compout_3", sep=" ")
    rr = pd.DataFrame()
    for i in range(0,len(df_3),21):
        rr = rr.append(pd.DataFrame(pd.concat([df_3.loc[i],df_3.loc[i+1],df_3.loc[i+2],df_3.loc[i+3],df_3.loc[i+4],df_3.loc[i+5],df_3.loc[i+6],df_3.loc[i+7],df_3.loc[i+8],df_3.loc[i+9],df_3.loc[i+10],df_3.loc[i+11],df_3.loc[i+12],df_3.loc[i+13],df_3.loc[i+14],df_3.loc[i+15],df_3.loc[i+16],df_3.loc[i+17],df_3.loc[i+18],df_3.loc[i+19],df_3.loc[i+20]],axis=0)).transpose()).reset_index(drop=True)
    rr.columns = head3
    cotrdi= pd.concat([zz,ss,rr],axis=1)
    cotrdi.to_csv(out, index=None, encoding='utf-8')
    os.remove('compout_1')
    os.remove('compout_2')
    os.remove('compout_3')
    #rr.to_csv(filename+".ctd_dist", index=None, encoding='utf-8')

def ctd_split(file,v,out):
    ctd(file,'tempfile_out')
    df4 = pd.read_csv("tempfile_out")
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('tempfile_out')

def ctd_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    ctd('input_sam.csv',result_filename)
    os.remove('input_sam.csv')
def ctd_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    ctd('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ctd_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    ctd('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ctd_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    ctd('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ctd_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    ctd('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def ctd_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    file1 =split('input_sam.csv',sp)
    file1.to_csv('sam_input.csv', index=None, header=False)
    ctd_split('sam_input.csv',sp,result_filename)
    os.remove('sam_input.csv')
    os.remove('input_sam.csv')
################################################################################
def searchAAIndex(AAIndex):
    found = -1;
    for i in range(len(AAIndexNames)):
        if(str(AAIndex) == AAIndexNames.iloc[i][0]):
            found = i;
    return found;
	
def phychem_AAI(file,AAIn,mode): 
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq = file;
    if(type(AAIn) == str):
        AAI = pd.read_csv(AAIn,header=None, sep=',');
        print (AAI.shape)
        AAI = AAI.values.tolist()
        print(AAI)
    else:
        AAI = AAIn;
    l2 = len(AAI)
    header  = AAI[0:l2];
    final=[];
    final.append(AAI);
    l1 = len(seq);
    seq=[seq[i].upper() for i in range(l1)];
    for i in range(l1):
        coded = encode(seq[i]);
        temp=[];
        for j in range(l2):
            pos = searchAAIndex(AAI[j]);
            sum=0;
            for k in range(len(coded)):
                val = AAIndex.iloc[pos,int(coded[k])]
                sum=sum+val;
            avg = round(sum/len(seq[i]),3);
            temp.append(avg);
        final.append(temp);
    if mode == 'all' :        
        file = open('AAIndex_all','w')        
    if mode == 'NT' :        
        file = open('AAIndex_NT','w')
    if mode == 'CT' :        
        file = open('AAIndex_CT','w')	
    if mode == 'rest' :        
        file = open('AAIndex_rest','w')		
    with file:
        writer = csv.writer(file);
        writer.writerows(final);
    return final;
def AAIndex_Phychem(filename,mode='all',m=0,n=0):
    time.sleep(1.5)
    seq = [];    
    if(type(filename) == str):
        seq1 = pd.read_csv(filename,header=None, sep=',');
        seq1 = pd.DataFrame(seq1[0].str.upper());
        [seq.append(seq1.iloc[i][0]) for i in range(len(seq1))];
    else:
        seq  = filename;
    AAIn = AAindices;
    if(type(AAIn) == str):
        AAI = pd.read_csv(AAIn,header=None, sep=',');
        l2 = AAI.shape[1]
        AAI = AAI.values.tolist()
        AAI  = list(itertools.chain.from_iterable(AAI))
    else:
        AAI = AAIn;
    l2 = len(AAI)
    l=len(seq)
    newseq=[];
    for i in range(0,l):
            l = len(seq[i]);
            if(mode=='NT'):
                n=m;
                print('Inside NT, m=',m,'n=',n)
                if(n!=0):
                    new = seq[i][0:n];
                    newseq.append(new);
                elif(n>l):
                    print('Warning! Sequence',i,"'s size is less than n. The output table would have NaN for this sequence");
                else:
                    print('Value of n is mandatory, it cannot be 0')
                    break;
            elif(mode=='CT'):
                n=m;
                if(n!=0):
                    new = seq[i][(len(seq[i])-n):]
                    newseq.append(new);
                elif(n>l):
                    print('WARNING: Sequence',i+1,"'s size is less than the value of n given. The output table would have NaN for this sequence");
                else:
                    print('Value of n is mandatory, it cannot be 0')
                    break;
            elif(mode=='all'):
                newseq = seq;
            elif(mode=='rest'):
                if(m==0):
                    print('Kindly provide start index for rest, it cannot be 0');
                    break;
                else:
                    #if(n<=len(seq[i])):
                    new = seq[i][m:len(seq[i])-n]
                    newseq.append(new)
            else:
                print("Wrong Mode. Enter 'NT', 'CT','all' or 'rest'");
    if(mode=='split'):
        newseq  = list(itertools.chain.from_iterable(newseq));
    phychem_AAI(newseq,AAI,mode);

def aai_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    AAIndex_Phychem('input_sam.csv','all')
    df = pd.read_csv('AAIndex_all')
    df.columns = 'AAI_'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('AAIndex_all')
    os.remove('input_sam.csv')
def aai_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    AAIndex_Phychem('input_sam.csv','NT',0,n)
    df = pd.read_csv('AAIndex_NT')
    df.columns = 'NAAI_'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('AAIndex_NT')
    os.remove('input_sam.csv')
def aai_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    AAIndex_Phychem('input_sam.csv','CT',0,c)
    df = pd.read_csv('AAIndex_CT')
    df.columns = 'CAAI_'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('AAIndex_CT')
    os.remove('input_sam.csv')
def aai_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    AAIndex_Phychem('input_sam.csv','rest',n,c)
    df = pd.read_csv('AAIndex_rest')
    df.columns = 'RAAI_'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('AAIndex_rest')
    os.remove('input_sam.csv')
def aai_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    AAIndex_Phychem('sam_input.csv','all')
    df = pd.read_csv('AAIndex_all')
    df.columns = 'NCAAI_'+df.columns
    df.to_csv(result_filename,index=None)
    os.remove('sam_input.csv')
    os.remove('AAIndex_all')
    os.remove('input_sam.csv')
def aai_st(file,out,v):
    readseq(file,'input_sam.csv')
    file1 = split('input_sam.csv',v)
    file2 = file1.iloc[:,0]
    AAIndex_Phychem(file2,'all')
    df4 = pd.read_csv('AAIndex_all')
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5.columns = 'AAI_'+df5.columns
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove('AAIndex_all')
    os.remove('input_sam.csv')


##################Parts of Seqeunces#########################
def nt(file,n):
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        df3.append(df2[0][i][0:n])
        df4 = pd.DataFrame(df3)
        #df4.to_csv(filename+".nt", index = None, header = False, encoding = 'utf-8')
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss < n:
            print('\nSequence number',i+1,'has length of',ss,'which is less than the provided value of N-terminal, that is',n,'. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
    return df4
def ct(file,n):
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        df3.append(df2[0][i][-n:])
        df4 = pd.DataFrame(df3)
        #df4.to_csv(filename+".ct", index = None, header = False, encoding = 'utf-8')
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss < n:
            print('\nSequence number',i+1,'has length of',ss,'which is less than the provided value of C-terminal, that is',n,'. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
    return df4
def rest(file,n,c):
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        if c == 0:
            df3.append(df2[0][i][n:])
        else:
            df3.append(df2[0][i][n:-c])
    df4 = pd.DataFrame(df3)
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss == 0:
            print('\nSequence number',i+1,'has length of',ss,'after removing provided N- and C-terminal residues, that is',str(n)+','+str(c),'. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
        #df4.to_csv(filename+".rest", index = None, header = False, encoding = 'utf-8')
    return df4
def restnc(file,n):
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        if n == 0:
            df3.append(df2[0][i][n:])
        else:
            df3.append(df2[0][i][n:-n])
    df4 = pd.DataFrame(df3)
    df4 = pd.DataFrame(df3)
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss == 0:
            print('\nSequence number',i+1,'has length of',ss,'after removing provided',str(n),' residues from N- and C-terminal. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
        #df4.to_csv(filename+".rest", index = None, header = False, encoding = 'utf-8')
    return df4
def nct(file,n):
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        df3.append(df2[0][i][:n]+df2[0][i][::-1][:n])
        df4 = pd.DataFrame(df3)
        #df4.to_csv(filename+".ct", index = None, header = False, encoding = 'utf-8')
    for i in range(0,len(df4)):
        ss = len(df4[0][i])
        if ss/2 < n:
            print('\nSequence number',i+1,'has length of',int(ss/2),'which is less than the provided value of N- and C-terminal, that is',n,'. Kindly check the sequences.')
            os.remove(file_output)
            sys.exit()
    return df4
def split(file,v):
    df1 = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df1[0].str.upper())
    for i in range(0,len(df2)):
        ss = len(df2[0][i])
        if ss == 0:
            print('\nSequence number',i+1,'has length of',ss,'. Hence, the number of splits should be between 2 to',ss,'. Kindly provide number of splits in the suggested range.')
            os.remove(file_output)
            sys.exit()
    k1 = []
    for e in range(0,len(df2)):
        s = 0
        k2 = []
        r = 0
        if len(df2[0][e])%v == 0:
            k2.extend(repeat(int(len(df2[0][e])/v),v))
        else:
            r = int(len(df2[0][e])%v)
            k2.extend(repeat(int(len(df2[0][e])/v),v-1))
            k2.append((int(len(df2[0][e])/v))+r)
        for j in k2:
            df3 = df2[0][e][s:j+s]
            k1.append(df3)
            s = j+s
    df4 = pd.DataFrame(k1)
    return df4
###########################################PSSM_Processing######################################
def pssm_comp(file,out):
   filename, file_ext = os.path.splitext(file)
   aa = list("ACDEFGHIKLMNPQRSTVWY")
   df=pd.read_csv(file, header=None)
   Matrix = [[0 for x in range(0,20)] for y in range(0,20)]
   j = 0
   df2 = []
   for i in aa:
      if i in list(df[0]):
         df1 = df.loc[df[0]==i]
         df2 = (df1.iloc[:,1:].sum(axis=0)/len(df))
         Matrix[j] =np.asarray(df2)
         j = j+1
      else:
         j = j+1
   f = open(out, 'w')
   sys.stdout = f
   for each in Matrix:
       for j in each:
           print("%.4f"%j,end=",")
   f.close()

def pssm_n1(file,out):
   df=pd.read_csv(file,header=None)
   df1 = df.iloc[:,0:21]
   def pssm(x):
       if type(x) is str:
           return x
       elif x:
           return (1/(1+(2.7182)**(-x)))
       elif x == -32768:
            return 0
       else:
            return
   df2 = df1.applymap(pssm)
   df2.to_csv(out,encoding='utf-8',index=False, header=False)
def pssm_n2(file,out):
   df=pd.read_csv(file,header=None)
   df1 = df.iloc[:,0:21]
   def pssm(x):
       a = 1000
       b = -1000
       if type(x) is str:
           return x
       elif x:
           return (x-a)/(b - a)
       else:
           return
   df2 = df1.applymap(pssm)
   df2.to_csv(out,encoding='utf-8',index=False, header=False)
def pssm_n3(file,out):
   df=pd.read_csv(file,header=None)
   df1 = df.iloc[:,0:21]
   def pssm(x):
       a = 1000
       b = -1000
       if type(x) is str:
           return x
       elif x:
           return ((x-a)*100)/(b - a)
       else:
           return
   df2 = df1.applymap(pssm)
   df2.to_csv(out,encoding='utf-8',index=False, header=False)
def pssm_n4(file,out):
   df=pd.read_csv(file,header=None)
   df1 = df.iloc[:,0:21]
   def pssm(x):
       if type(x) is str:
           return x
       elif x:
           return (1/(1+(2.7182)**(-x/100)))
       else:
           return
   df2 = df1.applymap(pssm)
   df2.to_csv(out,encoding='utf-8',index=False, header=False)

def aab(file,out):
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    uu = []
    for ss in df[0]:
         uu.append(len(ss))
    zz = df.iloc[:,0]
    f = open(out, mode='w')
    sys.stdout = f
    A=('1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    C=('0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    D=('0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    E=('0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    F=('0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    G=('0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    H=('0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0')
    I=('0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0')
    K=('0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0')
    L=('0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0')
    M=('0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0')
    N=('0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0')
    P=('0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0')
    Q=('0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0')
    R=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0')
    S=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0')
    T=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0')
    V=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0')
    W=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0')
    Y=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1')
    for mm in range (1,max(uu)+1):
        for ee in std:
            print(ee+str(mm),end=',')
    print("")
    for i in range(0,len(zz)):
        for j in zz[i]:
            if j == "A":
                print(''.join(A), end = ',')
            if j == "C":
                print(''.join(C), end = ',')
            if j == "D":
                print(''.join(D), end = ',')
            if j == "E":
                print(''.join(E), end = ',')
            if j == "F":
                print(''.join(F), end = ',')
            if j == "G":
                print(''.join(G), end = ',')
            if j == "H":
                print(''.join(H), end = ',')
            if j == "I":
                print(''.join(I), end = ',')
            if j == "K":
                print(''.join(K), end = ',')
            if j == "L":
                print(''.join(L), end = ',')
            if j == "M":
                print(''.join(M), end = ',')
            if j == "N":
                print(''.join(N), end = ',')
            if j == "P":
                print(''.join(P), end = ',')
            if j == "Q":
                print(''.join(Q), end = ',')
            if j == "R":
                print(''.join(R), end = ',')
            if j == "S":
                print(''.join(S), end = ',')
            if j == "T":
                print(''.join(T), end = ',')
            if j == "V":
                print(''.join(V), end = ',')
            if j == "W":
                print(''.join(W), end = ',')
            if j == "Y":
                print(''.join(Y), end = ',')
        print("")
    f.truncate()
def aab_split(file,n,out):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    file1 = split(file,n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aab('sam_input.csv','tempfile_out')
    ff1 = open(out,'w')
    sys.stdout=ff1
    aa = []
    for i in file1[0]:
        aa.append(len(i))
    uu = []
    for i in range(0,len(aa),n):
        uu.append(aa[i:i+n])
    for i in range(1,n+1):
        kk = max(uu)[i-1]
        for j in range(1,kk+1):
            for k in std:
                print(k+str(j)+'_s'+str(i), end=",")
    print("")
    with open("tempfile_out","r") as f:
        fob = f.readlines()
        fob_1 = fob[1:]
    for each in range(0,len(fob_1),n):
        print(','.join(fob_1[each:each+n]).replace(",\n,",",").replace("\n",""))
    ff1.truncate()
    os.remove('sam_input.csv')
    os.remove('tempfile_out')

def bin_di(file,q,out):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    filename, file_extension = os.path.splitext(file)
    df2=pd.read_csv(file,header=None)
    df = pd.DataFrame(df2[0].str.upper())
    uu = []
    for ss in df[0]:
         uu.append(len(ss))
    mat3 = pd.read_csv("bin_di.csv", header = None)
    mat3.set_index(0, inplace = True)
    mat3.index = pd.Series(mat3.index).replace(np.nan,'NA')
    f = open(out, 'w')
    sys.stdout = f
    for ss in range (1,(max(uu)-q)+1):
        for uu in std:
            for mm in std:
                print(uu+mm+str(ss),end=",")
    print("")
    for i in range(0,len(df)):
        for j in range(0,(len(df[0][i])-(q))):
            temp1 = df[0][i][j:j+q+1:q]
            for each in (mat3.loc[temp1].values.ravel()):
                print("%.0f"%each, end = ",", flush = True)
        print("")
    f.truncate()

def dpb_split(file,n,q,out):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    file1 = split(file,n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bin_di('sam_input.csv',q,'tempfile_out')
    ff1 = open(out,'w')
    sys.stdout=ff1
    aa = []
    for i in file1[0]:
        aa.append(len(i))
    uu = []
    for i in range(0,len(aa),n):
        uu.append(aa[i:i+n])
    for i in range(1,n+1):
        kk = max(uu)[i-1]
        for j in range(1,(kk+1)-q):
            for k in std:
                for pp in std:
                    print(k+pp+str(j)+'_s'+str(i), end=",")
    print("")
    with open("tempfile_out","r") as f:
        fob = f.readlines()
        fob_1 = fob[1:]
    for each in range(0,len(fob_1),n):
        print(','.join(fob_1[each:each+n]).replace(",\n,",",").replace("\n",""))
    ff1.truncate()
    os.remove('sam_input.csv')
    os.remove('tempfile_out')

#########################################Atomic binary###############################
def atom_bin(file,out) :
    df=pd.read_csv(file,header=None)
    ############binary matrix for atoms
    f = open('matrix_atom.out', 'w')
    sys.stdout = f
    print("C,H,N,O,S,")
    x = []
    for i in range(0,5) :
        x.append([])
        for j in range(0,5) :
            if i == j :
                x[i].append(1)
            else :
                x[i].append(0)

            print(x[i][j], end=",")
        print("")
    f.truncate()
##############associate binary values to atoms
    mat = pd.read_csv("matrix_atom.out")
    mat1 = mat.iloc[:,:-1]
    mat2 = mat1.transpose()
    df1 = pd.read_csv(paths_2+"/atom.csv",header=None)
    zz = []
    kk = pd.DataFrame()
    df1 = pd.read_csv(paths_2+"/atom.csv",header=None)
    for i in range(0,len(df1)) :
        zz.append([])
        for j in range(0,len(df1[1][i])) :
            temp = df1[1][i][j]
            zz[i].append(mat2.loc[temp])

    f1 = open('bin_atom_1', 'w')
    sys.stdout = f1
    for i in range(0,len(zz)) :
        for row in zz[i]:
            print(",".join(map(str,row)), end=",")
        print("")

    f1.truncate()

    with open('bin_atom_1', 'r') as f:
        g = list(f)

    for i in range(0,len(g)) :
        g[i] = g[i].replace(",\n","")

    df1["bin"]=g

    #########binary atom values for given file
    ss =[]
    for i in range(0,len(df)):
        ss.append([])
        for j in df[0][i]:
            for k in range(0,len(df1)):
                if j==df1[0][k]:
                    ss[i].append(df1[1][k])
    uu = []
    for i in range (0,len(ss)):
        uu.append(len("".join(ss[i])))

    xx=[]
    jj = 0
    for i in range(0,len(df)) :
        xx.append([])
        while jj < len(df[0][i]) :
            temp=df[0][i][jj]
            for k in range(0,len(df1)) :
                if temp == df1[0][k][0] :
                    xx[i].append(df1.iloc[k,2])
            jj += 1
        jj = 0


    f2 = open(out, 'w')
    sys.stdout = f2
    for pp in range(1,max(uu)+1):
        for aa in ('C','H','N','O','S'):
            print(aa+str(pp),end=",")
    print("")

    for i in range(0,len(xx)) :
        for row in xx[i]:
            print("".join(map(str,row)), end=",")
        print("")
    f2.truncate()
    os.remove('bin_atom_1')
    os.remove('matrix_atom.out')

def atb_split(file,n,out):
    std_1 = list("CHNOS")
    file1 = split(file,n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atom_bin('sam_input.csv','tempfile_out')
    df123 = pd.read_csv('tempfile_out')
    df456 = df123.iloc[:,:-1]
    ff1 = open(out,'w')
    sys.stdout=ff1
    aa = []
    for i in range(0,len(df456)):
        aa.append(len(df456.loc[i])/5)
    uu = []
    for i in range(0,len(aa),n):
        uu.append(aa[i:i+n])
    for i in range(1,n+1):
        kk = max(uu)[i-1]
        for j in range(1,int(kk+1)):
            for k in std_1:
                print(k+str(j)+'_s'+str(i), end=",")
    print("")
    with open("tempfile_out","r") as f:
        fob = f.readlines()
        fob_1 = fob[1:]
    for each in range(0,len(fob_1),n):
        print(','.join(fob_1[each:each+n]).replace(",\n,",",").replace("\n",""))
    ff1.truncate()
    os.remove('sam_input.csv')
    os.remove('tempfile_out')

def bond_bin(file,out) :
    df=pd.read_csv(file,header=None)
    ############binary matrix for atoms
    f = open('matrix_can_pat.out', 'w')
    sys.stdout = f
    print("-,=,c,b,")
    x = []
    for i in range(0,4) :
        x.append([])
        for j in range(0,4) :
            if i == j :
                x[i].append(1)
            else :
                x[i].append(0)

            print(x[i][j], end=",")
        print("")
    f.truncate()

##############associate binary values to bonds
    mat = pd.read_csv("matrix_can_pat.out")
    mat1 = mat.iloc[:,:-1]
    mat2 = mat1.transpose()
    df1 = pd.read_csv(paths_2+"/can_pat.csv")
    zz = []
    kk = pd.DataFrame()

    for i in range(0,len(df1)) :
        zz.append([])
        for j in range(0,len(df1.iloc[:,1][i])) :
            temp = str(df1.iloc[:,1][i][j])
            zz[i].append(mat2.loc[temp])

    f1 = open('bin_bond_1', 'w')
    sys.stdout = f1
    for i in range(0,len(zz)) :
        for row in zz[i]:
            print(",".join(map(str,row)), end=",")
        print("")

    f1.truncate()

    with open('bin_bond_1', 'r') as f:
        g = list(f)

    for i in range(0,len(g)) :
        g[i] = g[i].replace(",\n","")

    df1["bin"] = g


    ss =[]
    for i in range(0,len(df)):
        ss.append([])
        for j in df[0][i]:
            for k in range(0,len(df1)):
                if j==df1['Name'][k]:
                    ss[i].append(df1['canonical_pattern'][k])
    uu = []
    for i in range (0,len(ss)):
        uu.append(len("".join(ss[i])))

    xx=[]
    jj = 0
    for i in range(0,len(df)) :
        xx.append([])
        while jj < len(df[0][i]) :
            temp=df[0][i][jj]
            for k in range(0,len(df1)) :
                if temp == df1.iloc[k,0][0] :
                    xx[i].append(df1.iloc[k,2])
            jj += 1
        jj = 0

    f2 = open(out, 'w')
    sys.stdout = f2
    for pp in range(1,max(uu)+1):
        for aa in ('SI','DO','CY','BE'):
            print(aa+str(pp),end=",")
    print("")
    for i in range(0,len(xx)) :
        for row in xx[i]:
            print("".join(map(str,row)), end=",")
        print("")
    f2.truncate()
    os.remove('matrix_can_pat.out')
    os.remove('bin_bond_1')

def btb_split(file,n,out):
    std_1 = ['SI','DO','CY','BE']
    file1 = split(file,n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bond_bin('sam_input.csv','tempfile_out')
    df123 = pd.read_csv('tempfile_out')
    df456 = df123.iloc[:,:-1]
    ff1 = open(out,'w')
    sys.stdout=ff1
    aa = []
    for i in range(0,len(df456)):
        aa.append(len(df456.loc[i])/4)
    uu = []
    for i in range(0,len(aa),n):
        uu.append(aa[i:i+n])
    for i in range(1,n+1):
        kk = max(uu)[i-1]
        for j in range(1,int(kk+1)):
            for k in std_1:
                print(k+str(j)+'_s'+str(i), end=",")
    print("")
    with open("tempfile_out","r") as f:
        fob = f.readlines()
        fob_1 = fob[1:]
    for each in range(0,len(fob_1),n):
        print(','.join(fob_1[each:each+n]).replace(",\n,",",").replace("\n",""))
    ff1.truncate()
    os.remove('sam_input.csv')
    os.remove('tempfile_out')

def pcp_bin(file,out):
    df = pd.read_csv(file, header = None)
    uu = []
    for ss in df[0]:
         uu.append(len(ss))
    zz = df.iloc[:,0]
    f = open(out, mode='w')
    sys.stdout = f
    A=('0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,0,1,0,0,1,1,0')
    C=('0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,1,0,0,1,1,0')
    D=('0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0')
    E=('0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1')
    F=('0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,1')
    G=('0,0,1,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,1,0')
    H=('1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1')
    I=('0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,1')
    K=('1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1')
    L=('0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,1')
    M=('0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,1,0,0,1')
    N=('0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,1,0')
    P=('0,0,1,0,1,1,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,1,0,1,0')
    Q=('0,0,1,1,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,1')
    R=('1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1')
    S=('0,0,1,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,1,1,1,0')
    T=('0,0,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,0,0,1,0,1,0')
    V=('0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,1,0')
    W=('0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,1')
    Y=('0,0,1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,1')
    for mm in range (1,max(uu)+1):
        for ee in ('PC','NC','NE','PO','NP','AL','CY','AR','AC','BS','NE_pH','HB','HL','NT','HX','SC','SS_HE','SS_ST','SS_CO','SA_BU','SA_EX','SA_IN','TN','SM','LR'):
            print(ee+str(mm),end=',')
    print("")
    for i in range(0,len(zz)):
        for j in zz[i]:
            if j == "A":
                print(''.join(A), end = ',')
            if j == "C":
                print(''.join(C), end = ',')
            if j == "D":
                print(''.join(D), end = ',')
            if j == "E":
                print(''.join(E), end = ',')
            if j == "F":
                print(''.join(F), end = ',')
            if j == "G":
                print(''.join(G), end = ',')
            if j == "H":
                print(''.join(H), end = ',')
            if j == "I":
                print(''.join(I), end = ',')
            if j == "K":
                print(''.join(K), end = ',')
            if j == "L":
                print(''.join(L), end = ',')
            if j == "M":
                print(''.join(M), end = ',')
            if j == "N":
                print(''.join(N), end = ',')
            if j == "P":
                print(''.join(P), end = ',')
            if j == "Q":
                print(''.join(Q), end = ',')
            if j == "R":
                print(''.join(R), end = ',')
            if j == "S":
                print(''.join(S), end = ',')
            if j == "T":
                print(''.join(T), end = ',')
            if j == "V":
                print(''.join(V), end = ',')
            if j == "W":
                print(''.join(W), end = ',')
            if j == "Y":
                print(''.join(Y), end = ',')
        print("")
    f.truncate()
def pcb_split(file,n,out):
    std_3 = ['PC','NC','NE','PO','NP','AL','CY','AR','AC','BS','NE_pH','HB','HL','NT','HX','SC','SS_HE','SS_ST','SS_CO','SA_BU','SA_EX','SA_IN','TN','SM','LR']
    file1 = split(file,n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_bin('sam_input.csv','tempfile_out')
    ff1 = open(out,'w')
    sys.stdout=ff1
    aa = []
    for i in file1[0]:
        aa.append(len(i))
    uu = []
    for i in range(0,len(aa),n):
        uu.append(aa[i:i+n])
    for i in range(1,n+1):
        kk = max(uu)[i-1]
        for j in range(1,kk+1):
            for k in std_3:
                print(k+str(j)+'_s'+str(i), end=",")
    print("")
    with open("tempfile_out","r") as f:
        fob = f.readlines()
        fob_1 = fob[1:]
    for each in range(0,len(fob_1),n):
        print(','.join(fob_1[each:each+n]).replace(",\n,",",").replace("\n",""))
    ff1.truncate()
    os.remove('sam_input.csv')
    os.remove('tempfile_out')

def aai_bin(file,out):
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    uu = []
    for ss in df[0]:
         uu.append(len(ss))
    zz = df.iloc[:,0]
    f = open(out, mode='w')
    sys.stdout = f
    A=('0,0,1,1,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,1,1,1,1,1,0,1,0,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,1,0,1,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,0,1,0,1,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0,0,0,1,1,0,1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,1,1,1,1,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,1,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,1,0,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0')
    C=('1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,1,0,1,0,0,1,0,0,0,1,1,0,1,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,1,1,0,0,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,1,0,1,1,0,1,0,1,0,1,1,0,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,1,1,0,1,1,1,1,1,0,1,0,0,0,0,1,0,0,1,0,1,1,1,1,1,1,0,1,1,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,0,0,0,0,0,1,0,1,0,1,1,1,0,0,0,1,0,1,0,1,0,1,1,0,0,0,1,0,1,0,1,0,0,0,0,1,1,1,1,0,1,0,0,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,0')
    D=('1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,1,0,1,0,1,0,1,1,0,0,1,0,0,1,0,0,0,1,0,0,1,1,0,1,1,0,1,1,0,0,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,0,1,0,1,1,1,1,0,1,0,1,0,0,0,0,1,1,1,0,0,1,0,0,0,0,1,1,1,0,0,0,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,1,1,0,0,0,1,0,1,0,0,0,0,1,0,1,1,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,1,1,1,1,1,1,0,1,0,1,0,1,0,0,0,1,1,0,1,0,0,0,0,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,0,0,0,0,1,1,1,0,0,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,1,1,0,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,0,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,1,0,0')
    E=('0,0,0,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,1,0,0,1,1,0,1,1,1,0,1,0,0,0,1,1,1,1,0,1,1,1,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,1,1,0,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1,0,1,0,0,0,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,1,1,1,1,1,0,0,1,1,1,0,1,0,0,1,0,1,1,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,1,1,0,1,0,1,1,1,1,1,1,1,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,0,0,0,0,0,1,0,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,1,0,0,0,1,0,1,1,0,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,1,0,0')
    F=('1,1,1,1,1,1,0,0,1,0,0,1,1,0,1,0,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1,0,0,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,1,0,0,1,0,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,1,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,1,0,0,1,1,0,1,1,0,0,0,0,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,1,0,1,1,1,1,1,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,0,0,1,1,1,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,0,1,1,1,1,1,1,0,1,1,1,1,0,0,1,0,1,1,0,1,1,0,0,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,1,1,1,1,1,0,0,0,0,1,0,1,1,1,1,0,0,0,0,1,1,0,1,1,1,1,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,1,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,1,1,0,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1')
    G=('0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,1,1,1,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,1,1,1,0,0,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,1,1,0,1,1,0,1,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1,0,1,0,0,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,1,0,1,0,0,1,1,0,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,1,1,1,0,0,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0')
    H=('1,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,1,1,1,0,1,0,0,0,0,0,1,1,0,1,1,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,1,1,1,0,0,0,0,1,0,1,0,0,1,1,1,0,1,0,0,0,0,0,1,1,1,0,0,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,1,1,0,0,0,0,0,1,0,1,1,0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,0,0,0,0,1,0,1,1,0,1,0,1,1,0,1,0,1,1,0,0,0,1,1,1,1,0,1,1,1,1,1,0,0,1,1,0,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,0,0,1,0,1,1,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,1,0,1,0,1,1,0,0,0,0,1,1,0,1,1,1,0,0,1,0,1,0,0,0,1,1,0,0,0,1,1,0,0,0,0,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1,1,1,0,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,1,1,0,1,0,0,0,1,0,1,1,1,1,0,0,0,1,0,0,1,1,1,1,1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,1,1,1,1,1,1,0,0,0,1,1,0,0,1,0,1,1,1,0,0,1,1,1,1,0,1,1,1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,0')
    I=('0,1,1,1,1,0,0,0,0,1,0,1,1,1,1,1,1,1,0,0,0,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,1,0,1,1,1,0,1,1,0,1,0,0,1,1,1,1,0,1,1,1,1,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,1,1,1,1,1,0,0,0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,0,1,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,0,0,1,0,0,1,0,1,1,0,0,0,1,0,1,0,1,1,0,1,0,1,0,1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,1,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,0,0,1,1,0,1,1,1,1,1,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,0,0,0,0,1,0,1,1,0,1,1,1,1,1,1,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1')
    K=('0,1,0,0,1,1,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,1,1,1,0,0,1,0,0,0,1,1,1,0,1,1,0,0,1,1,0,1,1,1,1,0,0,0,0,1,1,1,0,1,1,1,1,0,0,1,1,0,1,0,1,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,1,0,0,1,0,1,1,0,1,0,1,1,0,0,1,1,1,0,1,0,0,1,0,1,0,1,1,0,0,1,1,0,1,1,1,1,0,0,0,0,1,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,1,0,1,0,1,0,0,1,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,1,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1,0,1,1,0,1,1,0,1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,0,1,0,1,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,1,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,0,0,0,1,1,0,1,1,0,1,0,1,0,0,0,1,1,1,0,1,0,1,1,0,1,1,1,1,0,1,0,0,0,1,1,1,1,1,0,1,0,1,0,0,0,0,0,1,0,1,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,1,1,0')
    L=('0,1,1,1,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,1,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1,1,0,1,1,1,1,1,0,1,0,0,1,0,1,1,0,1,1,1,1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,0,1,1,1,0,1,0,1,1,0,0,1,1,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,1,0,0,0,0,1,0,0,1,1,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,0,1,1,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1')
    M=('1,1,1,1,0,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1,1,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,0,1,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,1,1,0,1,0,1,1,0,0,1,1,0,0,0,0,1,0,0,1,1,0,1,1,0,0,0,0,0,0,1,0,1,0,1,1,1,1,0,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,1,0,1,1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,1,1,0,1,0,0,0,0,0,1,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1')
    N=('1,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,0,1,0,0,1,1,0,1,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,1,0,0,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,1,1,0,1,0,0,1,1,1,0,0,0,1,1,0,1,1,1,1,1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,1,1,1,0,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,0,0,1,0,1,1,1,0,0,1,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,1,1,1,0,1,1,0,0,0,1,1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,0,1,1,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0')
    P=('1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,1,0,1,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,1,1,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,1,0,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,1,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,0,1,1,0,1,1')
    Q=('0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,1,0,0,0,0,0,1,0,1,1,0,1,0,1,1,0,0,1,1,1,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,1,0,1,0,1,1,1,0,0,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,1,1,0,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,1,1,0,1,0,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,1,0,0,1,0,1,0,1,1,1,0,1,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,0,1,0,1,1,1,1,1,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,1,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0')
    R=('0,0,0,0,1,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1,1,1,1,0,1,0,1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,0,0,0,0,0,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,1,1,0,1,0,0,1,1,1,0,1,1,0,0,1,1,0,1,1,0,1,0,0,0,0,1,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,0,1,1,0,1,1,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,0,1,1,1,1,1,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,0,0,1,0,0,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,1,0,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,1,1,0,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,1,1,0,0,1,1,1,1,1,1,0,1,1,0,1,0,1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,0,0,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,0')
    S=('1,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,1,1,1,0,1,0,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,1,1,0,1,1,1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,1,0,1,1,0,0,0,0,1,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,0,1,0,1,1,1,1,0,1,1,1,1,1,0,1,0,1,0,0,0,0,1,0,1,1,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,0,1,1,0,0,1,0,1,1,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,0,0,1,1,0,1,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,1,1,1,0,0,1,1,0,1,0,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0')
    T=('0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,1,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,1,0,1,0,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,1,0,0,0,0,0,0,1,1,1,1,1,0,0,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,1,1,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,1,1,1,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,0,1,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0')
    V=('0,1,1,1,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,1,0,0,1,1,1,1,0,1,1,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,1,0,1,0,0,1,0,0,1,0,0,1,1,1,1,0,1,0,0,0,0,0,1,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1,0,1,0,1,0,0,0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,1,1,0,0,1,1,0,0,1,1,1,1,1,1,0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,0,1,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,0,0,1,1,0,0,0,1,1,1,1,0,1,1,1,1,1,0,0,1')
    W=('1,1,0,1,0,0,1,0,0,1,1,1,1,0,0,1,1,0,1,0,0,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,1,0,0,0,0,0,1,1,1,0,1,0,1,1,1,0,1,1,1,0,0,0,0,1,1,0,1,0,0,0,1,1,0,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,0,1,1,0,1,1,0,1,1,1,1,1,0,1,1,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,0,1,0,1,1,1,1,0,0,0,0,1,0,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,0,1,1,1,0,1,0,0,1,0,1,1,0,1,1,0,1,0,0,1,1,1,1,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,0,1,1,0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,1,0,1,1,0,1,1,0,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,0,1,0,0,1,0,1,0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0,1,1,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,0,1,0,1,1,1,1,0,1,0,0,1')
    Y=('1,1,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,1,0,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,0,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,0,0,1,0,0,1,1,0,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,1,0,1,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,1,1,0,1,0,1,1,0,1,1,0,1,0,1,0,1,1,0,0,1,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,1,0,0,0,1,1,0,1,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,0,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0,0,1,1,1,1,1,0,1,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,0,0,1,0,0,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,1,0,1,0,1,1,1,0,1,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,1,0,0,1,0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,0,0,1')
    for mm in range (1,max(uu)+1):
        for ee in ('ANDN920101','ARGP820101','ARGP820102','ARGP820103','AURR980101','AURR980102','AURR980103','AURR980104','AURR980105','AURR980106','AURR980107','AURR980108','AURR980109','AURR980110','AURR980111','AURR980112','AURR980113','AURR980114','AURR980115','AURR980116','AURR980117','AURR980118','AURR980119','AURR980120','BAEK050101','BASU050101','BASU050102','BASU050103','BEGF750101','BEGF750102','BEGF750103','BHAR880101','BIGC670101','BIOV880101','BIOV880102','BLAM930101','BLAS910101','BROC820101','BROC820102','BULH740101','BULH740102','BUNA790101','BUNA790102','BUNA790103','BURA740101','BURA740102','CASG920101','CEDJ970101','CEDJ970102','CEDJ970103','CEDJ970104','CEDJ970105','CHAM810101','CHAM820101','CHAM820102','CHAM830101','CHAM830102','CHAM830103','CHAM830104','CHAM830105','CHAM830106','CHAM830107','CHAM830108','CHOC750101','CHOC760101','CHOC760102','CHOC760103','CHOC760104','CHOP780101','CHOP780201','CHOP780202','CHOP780203','CHOP780204','CHOP780205','CHOP780206','CHOP780207','CHOP780208','CHOP780209','CHOP780210','CHOP780211','CHOP780212','CHOP780213','CHOP780214','CHOP780215','CHOP780216','CIDH920101','CIDH920102','CIDH920103','CIDH920104','CIDH920105','COHE430101','CORJ870101','CORJ870102','CORJ870103','CORJ870104','CORJ870105','CORJ870106','CORJ870107','CORJ870108','COSI940101','COWR900101','CRAJ730101','CRAJ730102','CRAJ730103','DAWD720101','DAYM780101','DAYM780201','DESM900101','DESM900102','DIGM050101','EISD840101','EISD860101','EISD860102','EISD860103','ENGD860101','FASG760101','FASG760102','FASG760103','FASG760104','FASG760105','FASG890101','FAUJ830101','FAUJ880101','FAUJ880102','FAUJ880103','FAUJ880104','FAUJ880105','FAUJ880106','FAUJ880107','FAUJ880108','FAUJ880109','FAUJ880110','FAUJ880111','FAUJ880112','FAUJ880113','FINA770101','FINA910101','FINA910102','FINA910103','FINA910104','FODM020101','FUKS010101','FUKS010102','FUKS010103','FUKS010104','FUKS010105','FUKS010106','FUKS010107','FUKS010108','FUKS010109','FUKS010110','FUKS010111','FUKS010112','GARJ730101','GEIM800101','GEIM800102','GEIM800103','GEIM800104','GEIM800105','GEIM800106','GEIM800107','GEIM800108','GEIM800109','GEIM800110','GEIM800111','GEOR030101','GEOR030102','GEOR030103','GEOR030104','GEOR030105','GEOR030106','GEOR030107','GEOR030108','GEOR030109','GOLD730101','GOLD730102','GRAR740101','GRAR740102','GRAR740103','GUOD860101','GUYH850101','GUYH850102','GUYH850104','GUYH850105','HARY940101','HOPA770101','HOPT810101','HUTJ700101','HUTJ700102','HUTJ700103','ISOY800101','ISOY800102','ISOY800103','ISOY800104','ISOY800105','ISOY800106','ISOY800107','ISOY800108','JACR890101','JANJ780101','JANJ780102','JANJ780103','JANJ790101','JANJ790102','JOND750101','JOND750102','JOND920101','JOND920102','JUKT750101','JUNJ780101','JURD980101','KANM800101','KANM800102','KANM800103','KANM800104','KARP850101','KARP850102','KARP850103','KARS160101','KARS160102','KARS160103','KARS160104','KARS160105','KARS160106','KARS160107','KARS160108','KARS160109','KARS160110','KARS160111','KARS160112','KARS160113','KARS160114','KARS160115','KARS160116','KARS160117','KARS160118','KARS160119','KARS160120','KARS160121','KARS160122','KHAG800101','KIDA850101','KIMC930101','KLEP840101','KOEP990101','KOEP990102','KRIW710101','KRIW790101','KRIW790102','KRIW790103','KUHL950101','KUMS000101','KUMS000102','KUMS000103','KUMS000104','KYTJ820101','LAWE840101','LEVM760101','LEVM760102','LEVM760103','LEVM760104','LEVM760105','LEVM760106','LEVM760107','LEVM780101','LEVM780102','LEVM780103','LEVM780104','LEVM780105','LEVM780106','LEWP710101','LIFS790101','LIFS790102','LIFS790103','MANP780101','MAXF760101','MAXF760102','MAXF760103','MAXF760104','MAXF760105','MAXF760106','MCMT640101','MEEJ800101','MEEJ800102','MEEJ810101','MEEJ810102','MEIH800101','MEIH800102','MEIH800103','MITS020101','MIYS850101','MIYS990101','MIYS990102','MIYS990103','MIYS990104','MIYS990105','MONM990101','MONM990201','MUNV940101','MUNV940102','MUNV940103','MUNV940104','MUNV940105','NADH010101','NADH010102','NADH010103','NADH010104','NADH010105','NADH010106','NADH010107','NAGK730101','NAGK730102','NAGK730103','NAKH900101','NAKH900102','NAKH900103','NAKH900104','NAKH900105','NAKH900106','NAKH900107','NAKH900108','NAKH900109','NAKH900110','NAKH900111','NAKH900112','NAKH900113','NAKH920101','NAKH920102','NAKH920103','NAKH920104','NAKH920105','NAKH920106','NAKH920107','NAKH920108','NISK800101','NISK860101','NOZY710101','OLSK800101','ONEK900101','ONEK900102','OOBM770101','OOBM770102','OOBM770103','OOBM770104','OOBM770105','OOBM850101','OOBM850102','OOBM850103','OOBM850104','OOBM850105','PALJ810101','PALJ810102','PALJ810103','PALJ810104','PALJ810105','PALJ810106','PALJ810107','PALJ810108','PALJ810109','PALJ810110','PALJ810111','PALJ810112','PALJ810113','PALJ810114','PALJ810115','PALJ810116','PARJ860101','PARS000101','PARS000102','PLIV810101','PONJ960101','PONP800101','PONP800102','PONP800103','PONP800104','PONP800105','PONP800106','PONP800107','PONP800108','PONP930101','PRAM820101','PRAM820102','PRAM820103','PRAM900101','PRAM900102','PRAM900103','PRAM900104','PTIO830101','PTIO830102','PUNT030101','PUNT030102','QIAN880101','QIAN880102','QIAN880103','QIAN880104','QIAN880105','QIAN880106','QIAN880107','QIAN880108','QIAN880109','QIAN880110','QIAN880111','QIAN880112','QIAN880113','QIAN880114','QIAN880115','QIAN880116','QIAN880117','QIAN880118','QIAN880119','QIAN880120','QIAN880121','QIAN880122','QIAN880123','QIAN880124','QIAN880125','QIAN880126','QIAN880127','QIAN880128','QIAN880129','QIAN880130','QIAN880131','QIAN880132','QIAN880133','QIAN880134','QIAN880135','QIAN880136','QIAN880137','QIAN880138','QIAN880139','RACS770101','RACS770102','RACS770103','RACS820101','RACS820102','RACS820103','RACS820104','RACS820105','RACS820106','RACS820107','RACS820108','RACS820109','RACS820110','RACS820111','RACS820112','RACS820113','RACS820114','RADA880101','RADA880102','RADA880103','RADA880104','RADA880105','RADA880106','RADA880107','RADA880108','RICJ880101','RICJ880102','RICJ880103','RICJ880104','RICJ880105','RICJ880106','RICJ880107','RICJ880108','RICJ880109','RICJ880110','RICJ880111','RICJ880112','RICJ880113','RICJ880114','RICJ880115','RICJ880116','RICJ880117','ROBB760101','ROBB760102','ROBB760103','ROBB760104','ROBB760105','ROBB760106','ROBB760107','ROBB760108','ROBB760109','ROBB760110','ROBB760111','ROBB760112','ROBB760113','ROBB790101','ROSG850101','ROSG850102','ROSM880101','ROSM880102','ROSM880103','SIMZ760101','SNEP660101','SNEP660102','SNEP660103','SNEP660104','SUEM840101','SUEM840102','SUYM030101','SWER830101','TAKK010101','TANS770101','TANS770102','TANS770103','TANS770104','TANS770105','TANS770106','TANS770107','TANS770108','TANS770109','TANS770110','TSAJ990101','TSAJ990102','VASM830101','VASM830102','VASM830103','VELV850101','VENT840101','VHEG790101','VINM940101','VINM940102','VINM940103','VINM940104','WARP780101','WEBA780101','WERD780101','WERD780102','WERD780103','WERD780104','WILM950101','WILM950102','WILM950103','WILM950104','WIMW960101','WOEC730101','WOLR790101','WOLR810101','WOLS870101','WOLS870102','WOLS870103','YUTK870101','YUTK870102','YUTK870103','YUTK870104','ZASB820101','ZHOH040101','ZHOH040102','ZHOH040103','ZIMJ680101','ZIMJ680102','ZIMJ680103','ZIMJ680104','ZIMJ680105'):
            print(ee+'_'+str(mm),end=',')
    print("")
    for i in range(0,len(zz)):
        for j in zz[i]:
            if j == "A":
                print(''.join(A), end = ',')
            if j == "C":
                print(''.join(C), end = ',')
            if j == "D":
                print(''.join(D), end = ',')
            if j == "E":
                print(''.join(E), end = ',')
            if j == "F":
                print(''.join(F), end = ',')
            if j == "G":
                print(''.join(G), end = ',')
            if j == "H":
                print(''.join(H), end = ',')
            if j == "I":
                print(''.join(I), end = ',')
            if j == "K":
                print(''.join(K), end = ',')
            if j == "L":
                print(''.join(L), end = ',')
            if j == "M":
                print(''.join(M), end = ',')
            if j == "N":
                print(''.join(N), end = ',')
            if j == "P":
                print(''.join(P), end = ',')
            if j == "Q":
                print(''.join(Q), end = ',')
            if j == "R":
                print(''.join(R), end = ',')
            if j == "S":
                print(''.join(S), end = ',')
            if j == "T":
                print(''.join(T), end = ',')
            if j == "V":
                print(''.join(V), end = ',')
            if j == "W":
                print(''.join(W), end = ',')
            if j == "Y":
                print(''.join(Y), end = ',')
        print("")
    f.truncate()
def aib_split(file,n,out):
    std_3 = ['ANDN920101','ARGP820101','ARGP820102','ARGP820103','AURR980101','AURR980102','AURR980103','AURR980104','AURR980105','AURR980106','AURR980107','AURR980108','AURR980109','AURR980110','AURR980111','AURR980112','AURR980113','AURR980114','AURR980115','AURR980116','AURR980117','AURR980118','AURR980119','AURR980120','BAEK050101','BASU050101','BASU050102','BASU050103','BEGF750101','BEGF750102','BEGF750103','BHAR880101','BIGC670101','BIOV880101','BIOV880102','BLAM930101','BLAS910101','BROC820101','BROC820102','BULH740101','BULH740102','BUNA790101','BUNA790102','BUNA790103','BURA740101','BURA740102','CASG920101','CEDJ970101','CEDJ970102','CEDJ970103','CEDJ970104','CEDJ970105','CHAM810101','CHAM820101','CHAM820102','CHAM830101','CHAM830102','CHAM830103','CHAM830104','CHAM830105','CHAM830106','CHAM830107','CHAM830108','CHOC750101','CHOC760101','CHOC760102','CHOC760103','CHOC760104','CHOP780101','CHOP780201','CHOP780202','CHOP780203','CHOP780204','CHOP780205','CHOP780206','CHOP780207','CHOP780208','CHOP780209','CHOP780210','CHOP780211','CHOP780212','CHOP780213','CHOP780214','CHOP780215','CHOP780216','CIDH920101','CIDH920102','CIDH920103','CIDH920104','CIDH920105','COHE430101','CORJ870101','CORJ870102','CORJ870103','CORJ870104','CORJ870105','CORJ870106','CORJ870107','CORJ870108','COSI940101','COWR900101','CRAJ730101','CRAJ730102','CRAJ730103','DAWD720101','DAYM780101','DAYM780201','DESM900101','DESM900102','DIGM050101','EISD840101','EISD860101','EISD860102','EISD860103','ENGD860101','FASG760101','FASG760102','FASG760103','FASG760104','FASG760105','FASG890101','FAUJ830101','FAUJ880101','FAUJ880102','FAUJ880103','FAUJ880104','FAUJ880105','FAUJ880106','FAUJ880107','FAUJ880108','FAUJ880109','FAUJ880110','FAUJ880111','FAUJ880112','FAUJ880113','FINA770101','FINA910101','FINA910102','FINA910103','FINA910104','FODM020101','FUKS010101','FUKS010102','FUKS010103','FUKS010104','FUKS010105','FUKS010106','FUKS010107','FUKS010108','FUKS010109','FUKS010110','FUKS010111','FUKS010112','GARJ730101','GEIM800101','GEIM800102','GEIM800103','GEIM800104','GEIM800105','GEIM800106','GEIM800107','GEIM800108','GEIM800109','GEIM800110','GEIM800111','GEOR030101','GEOR030102','GEOR030103','GEOR030104','GEOR030105','GEOR030106','GEOR030107','GEOR030108','GEOR030109','GOLD730101','GOLD730102','GRAR740101','GRAR740102','GRAR740103','GUOD860101','GUYH850101','GUYH850102','GUYH850104','GUYH850105','HARY940101','HOPA770101','HOPT810101','HUTJ700101','HUTJ700102','HUTJ700103','ISOY800101','ISOY800102','ISOY800103','ISOY800104','ISOY800105','ISOY800106','ISOY800107','ISOY800108','JACR890101','JANJ780101','JANJ780102','JANJ780103','JANJ790101','JANJ790102','JOND750101','JOND750102','JOND920101','JOND920102','JUKT750101','JUNJ780101','JURD980101','KANM800101','KANM800102','KANM800103','KANM800104','KARP850101','KARP850102','KARP850103','KARS160101','KARS160102','KARS160103','KARS160104','KARS160105','KARS160106','KARS160107','KARS160108','KARS160109','KARS160110','KARS160111','KARS160112','KARS160113','KARS160114','KARS160115','KARS160116','KARS160117','KARS160118','KARS160119','KARS160120','KARS160121','KARS160122','KHAG800101','KIDA850101','KIMC930101','KLEP840101','KOEP990101','KOEP990102','KRIW710101','KRIW790101','KRIW790102','KRIW790103','KUHL950101','KUMS000101','KUMS000102','KUMS000103','KUMS000104','KYTJ820101','LAWE840101','LEVM760101','LEVM760102','LEVM760103','LEVM760104','LEVM760105','LEVM760106','LEVM760107','LEVM780101','LEVM780102','LEVM780103','LEVM780104','LEVM780105','LEVM780106','LEWP710101','LIFS790101','LIFS790102','LIFS790103','MANP780101','MAXF760101','MAXF760102','MAXF760103','MAXF760104','MAXF760105','MAXF760106','MCMT640101','MEEJ800101','MEEJ800102','MEEJ810101','MEEJ810102','MEIH800101','MEIH800102','MEIH800103','MITS020101','MIYS850101','MIYS990101','MIYS990102','MIYS990103','MIYS990104','MIYS990105','MONM990101','MONM990201','MUNV940101','MUNV940102','MUNV940103','MUNV940104','MUNV940105','NADH010101','NADH010102','NADH010103','NADH010104','NADH010105','NADH010106','NADH010107','NAGK730101','NAGK730102','NAGK730103','NAKH900101','NAKH900102','NAKH900103','NAKH900104','NAKH900105','NAKH900106','NAKH900107','NAKH900108','NAKH900109','NAKH900110','NAKH900111','NAKH900112','NAKH900113','NAKH920101','NAKH920102','NAKH920103','NAKH920104','NAKH920105','NAKH920106','NAKH920107','NAKH920108','NISK800101','NISK860101','NOZY710101','OLSK800101','ONEK900101','ONEK900102','OOBM770101','OOBM770102','OOBM770103','OOBM770104','OOBM770105','OOBM850101','OOBM850102','OOBM850103','OOBM850104','OOBM850105','PALJ810101','PALJ810102','PALJ810103','PALJ810104','PALJ810105','PALJ810106','PALJ810107','PALJ810108','PALJ810109','PALJ810110','PALJ810111','PALJ810112','PALJ810113','PALJ810114','PALJ810115','PALJ810116','PARJ860101','PARS000101','PARS000102','PLIV810101','PONJ960101','PONP800101','PONP800102','PONP800103','PONP800104','PONP800105','PONP800106','PONP800107','PONP800108','PONP930101','PRAM820101','PRAM820102','PRAM820103','PRAM900101','PRAM900102','PRAM900103','PRAM900104','PTIO830101','PTIO830102','PUNT030101','PUNT030102','QIAN880101','QIAN880102','QIAN880103','QIAN880104','QIAN880105','QIAN880106','QIAN880107','QIAN880108','QIAN880109','QIAN880110','QIAN880111','QIAN880112','QIAN880113','QIAN880114','QIAN880115','QIAN880116','QIAN880117','QIAN880118','QIAN880119','QIAN880120','QIAN880121','QIAN880122','QIAN880123','QIAN880124','QIAN880125','QIAN880126','QIAN880127','QIAN880128','QIAN880129','QIAN880130','QIAN880131','QIAN880132','QIAN880133','QIAN880134','QIAN880135','QIAN880136','QIAN880137','QIAN880138','QIAN880139','RACS770101','RACS770102','RACS770103','RACS820101','RACS820102','RACS820103','RACS820104','RACS820105','RACS820106','RACS820107','RACS820108','RACS820109','RACS820110','RACS820111','RACS820112','RACS820113','RACS820114','RADA880101','RADA880102','RADA880103','RADA880104','RADA880105','RADA880106','RADA880107','RADA880108','RICJ880101','RICJ880102','RICJ880103','RICJ880104','RICJ880105','RICJ880106','RICJ880107','RICJ880108','RICJ880109','RICJ880110','RICJ880111','RICJ880112','RICJ880113','RICJ880114','RICJ880115','RICJ880116','RICJ880117','ROBB760101','ROBB760102','ROBB760103','ROBB760104','ROBB760105','ROBB760106','ROBB760107','ROBB760108','ROBB760109','ROBB760110','ROBB760111','ROBB760112','ROBB760113','ROBB790101','ROSG850101','ROSG850102','ROSM880101','ROSM880102','ROSM880103','SIMZ760101','SNEP660101','SNEP660102','SNEP660103','SNEP660104','SUEM840101','SUEM840102','SUYM030101','SWER830101','TAKK010101','TANS770101','TANS770102','TANS770103','TANS770104','TANS770105','TANS770106','TANS770107','TANS770108','TANS770109','TANS770110','TSAJ990101','TSAJ990102','VASM830101','VASM830102','VASM830103','VELV850101','VENT840101','VHEG790101','VINM940101','VINM940102','VINM940103','VINM940104','WARP780101','WEBA780101','WERD780101','WERD780102','WERD780103','WERD780104','WILM950101','WILM950102','WILM950103','WILM950104','WIMW960101','WOEC730101','WOLR790101','WOLR810101','WOLS870101','WOLS870102','WOLS870103','YUTK870101','YUTK870102','YUTK870103','YUTK870104','ZASB820101','ZHOH040101','ZHOH040102','ZHOH040103','ZIMJ680101','ZIMJ680102','ZIMJ680103','ZIMJ680104','ZIMJ680105']
    file1 = split(file,n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aai_bin('sam_input.csv','tempfile_out')
    ff1 = open(out,'w')
    sys.stdout=ff1
    aa = []
    for i in file1[0]:
        aa.append(len(i))
    uu = []
    for i in range(0,len(aa),n):
        uu.append(aa[i:i+n])
    for i in range(1,n+1):
        kk = max(uu)[i-1]
        for j in range(1,kk+1):
            for k in std_3:
                print(k+str(j)+'_s'+str(i), end=",")
    print("")
    with open("tempfile_out","r") as f:
        fob = f.readlines()
        fob_1 = fob[1:]
    for each in range(0,len(fob_1),n):
        print(','.join(fob_1[each:each+n]).replace(",\n,",",").replace("\n",""))
    ff1.truncate()
    os.remove('sam_input.csv')
    os.remove('tempfile_out')


def aab_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    aab('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def aab_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aab('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def aab_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aab('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def aab_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aab('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def aab_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aab('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def aab_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    aab_split('input_sam.csv',sp,'temp_out')
    df = pd.read_csv('temp_out')
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    df2.to_csv('sam_allbin.aab_st',index=None)
    os.remove('temp_out')
    os.remove('input_sam.csv')

def dpb_wp(seq,result_filename,lg):
    readseq(seq,'input_sam.csv')
    bin_di('input_sam.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def dpb_nt(seq,result_filename,n,lg):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bin_di('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def dpb_ct(seq,result_filename,c,lg):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bin_di('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def dpb_rt(seq,result_filename,n,c,lg):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bin_di('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def dpb_nct(seq,result_filename,n,lg):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bin_di('sam_input.csv',lg,'tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def dpb_st(seq,result_filename,sp,lg):
    readseq(seq,'input_sam.csv')
    dpb_split('input_sam.csv',sp,lg,'temp_out')
    df = pd.read_csv('temp_out')
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('temp_out')
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
def atb_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    atom_bin('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('input_sam.csv')
    os.remove('tempfile_out')
def atb_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atom_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
def atb_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atom_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def atb_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atom_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def atb_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    atom_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def atb_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    atb_split('input_sam.csv',sp,'temp_out')
    df = pd.read_csv('temp_out')
    df2 = df.dropna(axis=1, how='all')
    df2 = df2.fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('temp_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def btb_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    bond_bin('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('input_sam.csv')
    os.remove('tempfile_out')
def btb_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bond_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
def btb_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bond_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def btb_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bond_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def btb_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    bond_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def btb_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    btb_split('input_sam.csv',sp,'temp_out')
    df = pd.read_csv('temp_out')
    df2 = df.dropna(axis=1, how='all')
    df2 = df2.fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('temp_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def pcb_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    pcp_bin('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('input_sam.csv')
    os.remove('tempfile_out')
def pcb_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
def pcb_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def pcb_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def pcb_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    pcp_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def pcb_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    pcb_split('input_sam.csv',sp,'temp_out')
    df = pd.read_csv('temp_out')
    df2 = df.dropna(axis=1, how='all')
    df2 = df2.fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('temp_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def aib_wp(seq,result_filename):
    readseq(seq,'input_sam.csv')
    aai_bin('input_sam.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('input_sam.csv')
    os.remove('tempfile_out')
def aib_nt(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nt('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aai_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'N'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
    os.remove('tempfile_out')
def aib_ct(seq,result_filename,c):
    readseq(seq,'input_sam.csv')
    file1 = ct('input_sam.csv',c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aai_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'C'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def aib_rt(seq,result_filename,n,c):
    readseq(seq,'input_sam.csv')
    file1 = rest('input_sam.csv',n,c)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aai_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'R'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def aib_nct(seq,result_filename,n):
    readseq(seq,'input_sam.csv')
    file1 = nct('input_sam.csv',n)
    file1.to_csv('sam_input.csv', index=None, header=False)
    aai_bin('sam_input.csv','tempfile_out')
    df = pd.read_csv('tempfile_out')
    df.columns = 'NC'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('tempfile_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
def aib_st(seq,result_filename,sp):
    readseq(seq,'input_sam.csv')
    aib_split('input_sam.csv',sp,'temp_out')
    df = pd.read_csv('temp_out')
    df2 = df.dropna(axis=1, how='all')
    df2 = df2.fillna('NA')
    df2.to_csv(result_filename,index=None)
    os.remove('temp_out')
    os.remove('input_sam.csv')
    os.remove('sam_input.csv')
##################################################################################
def pat1_str(file,w):
    if w%2==0:
        raise Exception('window size should be an odd number. The value of window size provided: {}'.format(w))
        sys.exit()
    df1 = pd.read_csv(file, header =None)
    ext = str('X')*int((w-1)/2)
    aa = []
    for i in range(0,len(df1)):
        new_str = ext+df1[0][i]+ext
        aa.append(new_str)
    df2 = pd.DataFrame(aa)
    bb = []
    for j in range(0,len(df2)):
        for k in range(0,len(df2[0][j])):
            cc = df2[0][j][k:k+w]
            if len(cc) == w:
                bb.append(cc)
#         df3 = pd.DataFrame(bb)
    return pd.DataFrame(bb)
def lookup_2(peptide,featureNum):
    PCP= pd.read_csv(paths_2+'/PhysicoChemical_X.csv', header=None) #Our reference table for properties
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = encode(peptide);

    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);


def pcp_2(file,outt):
    PCP= pd.read_csv(paths_2+'/PhysicoChemical_X.csv', header=None) #Our reference table for properties
    headers = ['PC','NC','NE','PO','NP','AL','CY','AR','AC','BS','NE_pH','HB','HL','NT','HX','SC','SS_HE','SS_ST','SS_CO','SA_BU','SA_EX','SA_IN','TN','SM','LR','Z1','Z2','Z3','Z4','Z5'];
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq  = file;

    l = len(seq);

    rows = PCP.shape[0]; # Number of features in our reference table
    col = 20 ; # Denotes the 20 amino acids

    seq=[seq[0][i].upper() for i in range(l)]
    sequenceFeature = [];
    sequenceFeature.append(headers); #To put property name in output csv

    for i in range(l): # Loop to iterate over each sequence
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): #Loop to iterate over each feature
            featureVal = lookup_2(seq[i],j)
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(round(featureVal/len(seq[i]),3));
            else:
                sequenceFeatureTemp.append('NaN')

        sequenceFeature.append(sequenceFeatureTemp);

    out = pd.DataFrame(sequenceFeature);

    file = open(outt,'w')
    with file:
        writer = csv.writer(file);
        writer.writerows(sequenceFeature);
    return sequenceFeature;
'''
Function Name: pat_pcp
Description: Gives 30 physico-chemical properties of a sequence

Function prototype: pat_pcp(file,w,outt)

Input:
@file: an input csv file with multiple sequences
@mode(optional, default = 'all'):
    Values possible:
        1) (default)'all' : to get features of entire protein
        2) 'NT' : to get the features of first n N-Terminal residues
        3) 'CT' : to get the features of last n C-Terminal residues
        4) 'rest' : to get the features of a sub-sequence from mth position to nth position(both inclusive)
@m(optional(mandatory if 'rest' is chosen, default=0): m is the start position of residue
@n(optional, default = '0'): n is the number of residues you want from desired terminal or end point (if 'rest' is chosen)


Output:
A matrix (csv file) of dimension (m x 30) containing sequences as rows and their 30 physico-chemical properties as columns
where m = number of sequences (each sequence separated by comma)

'''


'''

Function Name: pat_pcp
Description: Gives 30 physico-chemical properties of a sequence

Function prototype: pat_pcp(file,w,outt)

Input:
@file: an input csv file with multiple sequences
@mode(optional, default = 'all'):
    Values possible:
        1) (default)'all' : to get features of entire protein
        2) 'NT' : to get the features of first n N-Terminal residues
        3) 'CT' : to get the features of last n C-Terminal residues
        4) 'rest' : to get the features of a sub-sequence from mth position to nth position(both inclusive)
@m(optional(mandatory if 'rest' is chosen, default=0): m is the start position of residue
@n(optional, default = '0'): n is the number of residues you want from desired terminal or end point (if 'rest' is chosen)


Output:
A matrix (csv file) of dimension (m x 30) containing sequences as rows and their 30 physico-chemical properties as columns
where m = number of sequences (each sequence separated by comma)

'''


def pat_pcp(file,outt,w):
    m = 0
    n = 0
    mode = 'all'
    readseq(file,'input_sam.csv')
    file = pat1_str('input_sam.csv',w)
    if(type(file) == str):
        seq = pd.read_csv(file,header=None);
        #seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    #elif(type(file)==str && len(file))
    else:
        seq  = file;
    l = len(seq);

    newseq = [""]*l; # To store the trimmed sequence
    #print('Original Sequence');
    #print(seq)


    for i in range(0,l):

        #if(n<len(seq[i])):
        l = len(seq[0][i]);

        if(mode=='NT'):
            n=m;
            if(n!=0):
                newseq[i] = seq[i][0:n];

            elif(n>l):
                print('Warning! Sequence',i,"'s size is less than n. The output table would have NaN for this sequence");


            else:
                print('Value of n is mandatory, it cannot be 0')
                break;

        elif(mode=='CT'):
            n=m;
            if(n!=0):
                newseq[i] = seq[i][(len(seq[i])-n):]

            elif(n>l):
                print('WARNING: Sequence',i+1,"'s size is less than the value of n given. The output table would have NaN for this sequence");


            else:
                print('Value of n is mandatory, it cannot be 0')
                break;

        elif(mode=='all'):
            newseq = seq;

        elif(mode=='rest'):
            if(m==0):
                print('Kindly provide start index for rest, it cannot be 0');
                break;

            else:
                if(n<=len(seq[i])):
                    newseq[i] = seq[i][m:(len(seq[i])-n)]
                '''elif(n>len(seq[i])):
                    newseq[i] = seq[i][m-1:len(seq[i])]
                    print('WARNING: Since input value of n for sequence',i+1,'is greater than length of the protein, entire sequence starting from m has been considered')'''
        else:
            print("Wrong Mode. Enter 'NT', 'CT','all' or 'rest'");
    output = pcp_2(newseq,outt);
    return output
    os.remove('input_sam.csv')
def pattern(file,windowfile):
    with open(file,'r') as f:
      g = list(f)
    orig_stdout = sys.stdout
    n = open("input_sam.pat",'w')
    sys.stdout = n
    k= (int(windowfile)-1)/2
    new_str = "X" * int(k)
    for i in g:
          c = i.replace('\n','')
          c =new_str+c+new_str
          for j in range(0,len(c)):
              d = c[j:j+int(windowfile)]
              if len(d)==int(windowfile):
                  print(d.upper())
    n.truncate()
def pat_bin(file,output,windowfile):
    std_1 = list('ACDEFGHIKLMNPQRSTVWYX')
    readseq(file,'input_sam.csv')
    pattern('input_sam.csv',windowfile)
    df1 = pd.read_csv('input_sam.pat', header = None)
    df = pd.DataFrame(df1[0].str.upper())
    zz = df.iloc[:,0]
    f = open(output, 'w')
    sys.stdout = f
    A=('1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    C=('0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    D=('0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    E=('0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    F=('0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    G=('0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    H=('0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    I=('0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0')
    K=('0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0')
    L=('0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0')
    M=('0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0')
    N=('0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0')
    P=('0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0')
    Q=('0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0')
    R=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0')
    S=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0')
    T=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0')
    V=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0')
    W=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0')
    Y=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0')
    X=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    for e in range(1,windowfile+1):
        for t in std_1:
            print(t+str(e),end=",")
    print("")
    for i in range(0,len(zz)):
        for j in zz[i]:
            if j == "A":
                print(''.join(A), end = ',')
            if j == "C":
                print(''.join(C), end = ',')
            if j == "D":
                print(''.join(D), end = ',')
            if j == "E":
                print(''.join(E), end = ',')
            if j == "F":
                print(''.join(F), end = ',')
            if j == "G":
                print(''.join(G), end = ',')
            if j == "H":
                print(''.join(H), end = ',')
            if j == "I":
                print(''.join(I), end = ',')
            if j == "K":
                print(''.join(K), end = ',')
            if j == "L":
                print(''.join(L), end = ',')
            if j == "M":
                print(''.join(M), end = ',')
            if j == "N":
                print(''.join(N), end = ',')
            if j == "P":
                print(''.join(P), end = ',')
            if j == "Q":
                print(''.join(Q), end = ',')
            if j == "R":
                print(''.join(R), end = ',')
            if j == "S":
                print(''.join(S), end = ',')
            if j == "T":
                print(''.join(T), end = ',')
            if j == "V":
                print(''.join(V), end = ',')
            if j == "W":
                print(''.join(W), end = ',')
            if j == "Y":
                print(''.join(Y), end = ',')
            if j == "X":
                print(''.join(X), end = ',')
        print("")
    f.truncate()
    os.remove('input_sam.csv')
    os.remove('input_sam.pat')

def pat_csv(file,out,n):
    filename,file_ext = os.path.splitext(file)
    df3 = pd.read_csv(file, header=None)
    ss = []
    f = open(out, 'w')
    sys.stdout = f
    for i in range(0,len(df3)):
        ss.append([i])
        for j in range(0,(len(df3.loc[i])-n+1)):
            ss[i].append(df3.loc[i][j:j+n].values)
    for i in range(0,len(ss)) :
        for j in range(1,len(ss[i])) :
            for each in ss[i][j][0:len(ss[i][j])] :
                print(each, end=",")
        print("")
    f.truncate()

def pat_str(inputfile,outputfile,windowfile,extensionfile='y'):
    readseq(inputfile,'input_sam.csv')
    with open('input_sam.csv','r') as f:
        g = list(f)
    orig_stdout = sys.stdout
    n = open(outputfile,'w')
    sys.stdout = n
    aa =('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
    k= (int(windowfile)-1)/2
    new_str = "X" * int(k)
    for i in g:
        c = i.replace('\n','')
        if (extensionfile == 'y'):
            c =new_str+c+new_str
        for j in range(0,len(c)):
            d = c[j:j+int(windowfile)]
            if len(d)==int(windowfile):
                print(d.upper())
        print("")
    n.truncate()


def phychem_AAI_1(file,AAIn,outt):
    if(type(file) == str):
        seq = pd.read_csv(file,header=None);
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq = file;

    if(type(AAIn) == str):
        AAI = pd.read_csv(AAIn,header=None, sep=',');
        l2 = AAI.shape[1]
        AAI.values.tolist()
    else:
        AAI = AAIn;
    l2 = AAI.shape[1]
    header =['Sequence'];
    header.extend(i for i in AAIn.loc[0])
    final=[];
    final.append(header);
    l1 = len(seq);
    seq=[seq[0][i].upper() for i in range(l1)];
    for i in range(l1):

        coded = encode(seq[i]);
        temp=[];
        for j in range(l2):
            pos = searchAAIndex(AAI[j][0]);
            #print(pos)
            if(j==0):
                temp.append(seq[i]);
            sum=0;
            for k in range(len(coded)):
                val = AAIndex.iloc[pos,int(coded[k])]
                sum=sum+val;
            avg = round(sum/len(coded),3);
            temp.append(avg);
        final.append(temp);
    file = open(outt,'w')
    with file:
        writer = csv.writer(file);
        writer.writerows(final);
    return final;
def pat_aai(file,outt,w):
    AAIn = AAindices
    file = pat1_str(file,w)
    mode = 'all'
    m = 0
    n = 0
    if(type(file) == str):
        seq = pd.read_csv(file,header=None, sep=',');
        seq=seq.T
        seq[0].values.tolist()
        seq=seq[0];
    else:
        seq = file;

    if(type(AAIn) == str):
        AAI = pd.read_csv(AAIn,header=None, sep=',');
        l2 = AAI.shape[1]
        AAI.values.tolist()
    else:
        AAI = AAIn;
    l2 = AAI.shape[1]
    l=len(seq)
    newseq=[];
    for i in range(0,l):
            l = len(seq[0][i]);
            if(mode=='NT'):
                n=m;
                if(n!=0):
                    new = seq[i][0:n];
                    newseq.append(new);

                elif(n>l):
                    print('Warning! Sequence',i,"'s size is less than n. The output table would have NaN for this sequence");


                else:
                    print('Value of n is mandatory, it cannot be 0')
                    break;

            elif(mode=='CT'):
                n=m;
                if(n!=0):
                    new = seq[i][(len(seq[i])-n):]
                    newseq.append(new);
                elif(n>l):
                    print('WARNING: Sequence',i+1,"'s size is less than the value of n given. The output table would have NaN for this sequence");


                else:
                    print('Value of n is mandatory, it cannot be 0')
                    break;

            elif(mode=='all'):
                newseq = seq;

            elif(mode=='rest'):
                if(m==0):
                    print('Kindly provide start index for rest, it cannot be 0');
                    break;

                else:
                    new = seq[i][m:len(seq[i])-n]
                    newseq.append(new)
                    '''elif(n>len(seq[i])):
                        newseq[i] = seq[i][m-1:len(seq[i])]
                        print('WARNING: Since input value of n for sequence',i+1,'is greater than length of the protein, entire sequence starting from m has been considered')'''

            elif(mode=='split'):
                n=m;
                new = split(seq[i],m);
                newseq.append(new);
            else:
                print("Wrong Mode. Enter 'NT', 'CT','all' or 'rest'");
    if(mode=='split'):
        newseq  = list(itertools.chain.from_iterable(newseq));
    phychem_AAI_1(newseq,AAI,outt);
####################################File reading###################################
def readseq(file,out):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(array[1:]).upper())
        seqid.append(name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append("Seq_"+str(i))
    for i in seq:
        if 'B' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'J' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'O' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'U' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'Z' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'X' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
    df4 = pd.DataFrame(seq)
    df4.to_csv(out,index=None,header=False)
##############################################################################
