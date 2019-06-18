#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 16:45:02 2019

@author: Sana Habib 

This code is written for the data collected at Waimea Valley by the CMAIKI group at 
Univeristy of Hawaii at Manoa. This will load data files, clean data, present graphs of raw
data, filter data, compute and analyze distance matrices. 

An input required is a tsv file containing an abundance matrix with a header row and column. 
"""

from playdata import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



        
    

def main():
    #qualAnal()
    
    
#    df = addMD()
#    df = df.iloc[0:-4, :]
#    print(str(df.shape))
#    print(df.max().max())
#    pres, total = summarize(df)
#    t = 12311*254
#    print(str(sum(pres)/t))
    df = normMD()
#    df = df.iloc[0:-4, :]
#    print(str(df.shape))
#    print(df.max().max())
#    pres, total = summarize(df)
#    print(str(sum(pres)/t))
    l = df.shape[1]
#    vals = df[df!=0].stack()
#    vals = vals.sort_values()
#    print(vals[36166])
#    print(36166/len(vals))
#    print(vals[40687])
#    print(40687/len(vals))
#    print(vals[44755])
#    print(44755/len(vals))
#    plt.hist(vals, bins = 100, histtype = 'step')
    
    meta = df.iloc[-4:, 0:]
    meta = meta.transpose()
    meta = meta.loc[meta['Host'] == 'Plant']
    meta = meta.loc[meta['SampleType'] == 'LichenThallus']
    df = df.iloc[0:-4, :]
    for sample in meta.index:
        col = df[str(sample)]
        col = col.apply(pd.to_numeric)
        print(col.max())
        idmax = col.idxmax()
        print(idmax)
    
#    fig = plt.figure(0)
#    tri = np.zeros((l,l))
#    tri[np.triu_indices(l, 1)] = BCDM
#    tri = tri + tri.transpose()
#    plt.imshow(tri, cmap='hot', interpolation='nearest')
#    plt.colorbar()
#    plt.xlabel('samples')
#    plt.ylabel('samples')
#    plt.title('Bray Curtis Dissimilarity Heat Map')
#    txt = '1 corresponds to no similarity'
#    fig.text(.5, -.05, '1 corresponds to no similarity', ha='center')
#    fig.savefig('DMHeatmap.jpg')
  
#    df = df.transpose()
#    df.sort_values(1, inplace = True)
#    df = df.transpose()
##    DM = BrayCurtis(df)
#    BCDM = BrayCurtisDM(df.iloc[0:-4, :])
    
#    
#    DM = pDM(df.iloc[0:-4, :], p = 2 )
#    
#    
#    
#    A = adjacency(DM, .83)
#    D = degree(A, l)
#    L = np.diag(D) - CtoSymM(A, l)
#    
##    plt.imshow(L, cmap='hot', interpolation='nearest')
##    plt.colorbar()
##    plt.xlabel('samples')
##    plt.ylabel('samples')
##    plt.title('Laplacian Heat Map')
##    plt.savefig('LaplacianHeatMap2.jpg')
#    
#
#    [eigvals, eigvecs] = np.linalg.eigh(L, UPLO = 'U')
#    vals = abs(eigvals) > 10**-12
##    print(vals)
#    ind = np.where(vals)
#    ind = ind[0][0]
#    lambda2 = eigvals[ind]
##    print(eigvals)
#    print(lambda2)
#    f = eigvecs[:, ind]
#    z = range(len(f))
#    plt.scatter(f, z)
##    plt.xlabel('Real Number Line')
##    plt.ylabel('Eigenvector entry number')
##    plt.title('Entries of the second eigenvector')
##    plt.savefig('EigvectorEmpty.jpg')
#    
#    clusters = f>0
#    i = np.where(clusters == False)
#    print(len(i[0]))
#    i = i[0]
#    print(i)
#    print(CtoSymM(DM, l)[i])
#    data = df.iloc[0:-4, i]
#    data = data.apply(pd.to_numeric)
#    idmax = data.idxmax()
#    print(idmax)
#    print(data.max())
#    print(data.loc[48])
    
#    plt.xlabel('samples')
#    plt.ylabel('samples')
#    plt.title('2-norm Connectivity Heat Map')
#    fig.text(.5, -.05, '1 corresponds to no similarity', ha='center')
#    fig.savefig('2normconnectivitymp.jpg')
#    
#    np.savetxt('dm.txt',DM)
#    np.savetxt('connect.txt',C)
    
'''
recall: found alpha such that 80% of 
'''

    
            
if __name__ == "__main__":
    main()
    


