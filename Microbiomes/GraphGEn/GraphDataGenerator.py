#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 19:35:34 2019

@author: sanashahid
"""

import pandas as pd
import numpy as np
import sys
import math 
import matplotlib.pyplot as plt
from scipy import stats



def graphfile(df, A, e ):
    
    name = list(df)
    name = ['"'+str(n)+'"' for n in name]
    Group = list(df.loc['Host'])
    group = []
    for g in Group:
        if g == "Plant":
            group.append('0')
        elif g == "Animal":
            group.append('1')
        elif g == "Fungus":
            group.append('2')
        else:
            group.append('3')
    sources = []
    targets = []
    weights = []
    k = 0
    for i in range(254):
        for j in range(i+1, 254):
            a = A[k]
            if( a != 0):
                sources.append(str(name[i]))
                targets.append(str(name[j]))
                weights.append((a))
            k+=1
                
    filename = 'graphdata'+str(int(e*10))+'.json'
    file = open(filename, 'w')
    file.write('{\n"nodes": [\n')
    for i in range(len(name)):
        if (i == len(name)-1):
            file.write('{"id": '+name[i]+', "group": '+str(group[i])+ '} \n')
        else:
            file.write('{"id": '+name[i]+', "group": '+str(group[i])+ '}, \n')
    file.write('], \n"links":[\n')
    for i in range(len(weights)):
        if(i == len(weights)-1):
            file.write('{"source": '+sources[i]+', "target": '+targets[i]+', "value": '+str(weights[i])+'}\n')

        else:
            file.write('{"source": '+sources[i]+', "target": '+targets[i]+', "value": '+str(weights[i])+'},\n')
    file.write(']\n}')
    file.close()
    
def adjacency(DM, e): # adjacency matrix
    C = [];
    for x in DM:
        if x >= e: C.append(x) # keep connections that are similar in exp(-d^2) form 1 is similar, 0 is far
        else: C.append(0)
    return C



def pDM(fdata, p = 2): 
    '''
    input: pandas dataframe
    Compute the p-norm Distance matrix 
    return condensed upper triangular distance matrix
    '''      
    l = fdata.shape[1]
    dm = []
    for i in range(l):
        for j in range(i+1, l):
            dm.append(np.linalg.norm(fdata.iloc[:,i] - fdata.iloc[:,j], ord = p))
    dm = np.array(dm)
    #############################################################
    #These are new edits
    stat = np.percentile(dm[np.nonzero(dm)], 60.374)
    s  = stat
    DM = [math.exp(-(d**2)/(s**2)) for d in dm]
    print(stats.skew(DM))
    plt.figure()
    plt.hist(DM, bins = 20)
    plt.xlabel('Augmented Weights')
    plt.ylabel('Number of Samples')
    plt.title('Distribution of Augmented Weights')
    plt.savefig('AugWeights.jpg')
    plt.figure()
    plt.hist(dm, bins = 28)
    plt.xlabel('Weights')
    plt.ylabel('Number of Samples')
    plt.title('Distribution of  Weights')
    plt.savefig('Weights.jpg')
    
    
    return DM

def normMD():
    md = meta()
    data = imp()
    data = data.div(data.sum(axis = 0))
    df = pd.concat([data, md], sort = False, join = 'inner', axis = 0)
    return df

def meta(file = 'CMAIKI_Metadata.csv'):
    metadata = pd.read_csv(file, header = 0, index_col = 0)
    cols = [ 'SampleType', 'Habitat', 'Host', 'Trophic']
    md = metadata[ cols]
    names = md.index.values
    newnames = []
    for name in names:
        n =  '10' + str(name[6:10])
        newnames.append((n))
    pd.options.mode.chained_assignment = None
    md.loc[:,'id'] = (newnames)
    md = md.set_index( 'id' )
    md = md.transpose()
    return md

def imp(file = 'abundance_table_annotated_ID100.tsv'):
    '''
    input: string corresponding to file name
    import data contained in tsv file 
    return a pandas datafram containing all the data wihtout indexing columns
    '''
    data = pd.read_table(file, header = 0, index_col = 0)
    names = data.index.values
    newnames = []
    for name in names:
        try:
            newnames.append(int(name[4:]))
        except:
            newnames.append(name[4:])
    data['id'] = (newnames)
    data = data.set_index( 'id' )
    data = data.apply(pd.to_numeric)
    return data

def main():
   df = normMD()
   DM = pDM(df.iloc[0:-4, :], p = 2 )
   for e in [.1]:#, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1]:
        C = adjacency(DM, e)
        graphfile(df, C, e)
        
if __name__ == "__main__":
    main()