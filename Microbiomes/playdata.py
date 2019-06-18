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

import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.cluster.hierarchy import dendrogram, linkage
import pandas as pd
import numpy as np
import math




def imp(file = 'GivenFiles/abundance_table_annotated_ID100.tsv'):
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


def summarize(data):
    '''
    input: pandas dataframe
    Go through data and count presence, non present, and total abundance.
    return dataframe, list of zeros, presence counts, and total sums.
    '''
    d = data.replace(0, np.nan)
    pres = (d.count(axis = 1)).values
    total = data.sum(axis = 1).values
    return pres, total

def filt(data,present, limit = 200, op = 'l'):
    '''
    input: pandas dataframe, a vertical axis limit, and op = 'u' or 'l' for upper or lower filtering
    filter data based on presence count. Keep data that appears more than 5 times.
    return new data
    '''
    if op == 'l':present = present < limit
    else:present = present > limit
    fdata = data[present]
    return fdata.reset_index(drop = True)


def BrayCurtisDM(fdata):
    '''
    input: pandas dataframe
    Compute the Bray Curtis Distance matrix
    return condensed upper triangular bray curtis dissimilarity matrix
    '''
    l = fdata.shape[1]
    BCDM = []
    for i in range(l):
        for j in range(i+1, l):
            BCDM.append(distance.braycurtis(fdata.iloc[:,i], fdata.iloc[:,j]))
    BCDM = np.array(BCDM)
    return BCDM

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
#    print(stats.skew(DM))
#    plt.figure()
#    plt.hist(DM)
#    plt.figure()
#    plt.hist(dm)
    return DM

def dend(DM, name):
    '''
    input: condensed distance matrix, name of output file
    Source: https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.cluster.hierarchy.linkage.html
    This computes the linkage clusters in data 
    saves the dendogram image to file
    '''
    plt.figure()
    z = linkage(DM)
    dendrogram(z)
    plt.savefig(name+'.png')
    
def linreg(x,y, name1 = "", name2 = "", ylim = 100000):
    '''
    input: takes in two arrays
    computes a linear regression for the two arrays and plots the data with the linear regression
    output: returns linear regression model
    '''
    fn = np.polyfit(x, y, 1)
    f =  np.poly1d(fn)
    plt.figure()
    plt.plot(x,y, 'ko', x, f(x), 'r-')  
    plt.ylim(top = ylim)           ################
    plt.xlabel(name1)
    plt.ylabel(name2)
    plt.title(name1 + ' vs ' + name2)
    plt.savefig('Output/'+name1 + 'vs' + name2 + 'Scatter.png') 
    plt.close()
    return f

def rep(data, quant = ''):
    '''
    input: dataframe
    this summarizes the data frame, computes an abundance histogram, presence 
    histogram, and scatterplot of the two data sets. 
    output:it returns the prsence and abundance vectors
    '''
    zero, pres, total = summarize(data)
    n1, c1, n2, c2, n = Hist(total, quant + 'abundance', 39)
    print("In the %s histogram There are %d samples (%f)with %s count less than %d"%( 
          quant + 'abundance', n1,n1/n,'abundance', c1  ))
    print("In the %s histogram There are %d samples (%f) with %s count less than %d"%( 
          quant + 'abundance', n2,n2/n,'abundance', c2  ))
    n1, c1, n2, c2, n =Hist(pres, quant + 'Presence', 39)
    print("In the %s histogram There are %d samples (%f) with %s count less than %d"%( 
          quant + 'Presence', n1,n1/n,'Presence', c1  ))
    print("In the %s histogram There are %d samples(%f) with %s count less than %d"%( 
          quant + 'Presence', n2, n2/n,'Presence', c2  ))
    scat(pres, total, quant + 'presence' , quant + 'abundance')
    return pres, total

def viz(data):
    '''
    input: dataframe of abundance data
    this looks at the initial data, filters on presence, then based on a linear regression
    output: filtered data
    '''
    pres, total = rep(data, 'Original')
    data = filt(data, pres, limit = 100)
    pres, total = rep(data, 'filterenonPresence')
    f = linreg(pres, total, 'FilteredPresence' , 'FilteredAbundance')
    fardata = filt(data, abs(total - f(pres)), limit = 5000, op = 'g')
    closedata = filt(data, abs(total - f(pres)), limit = 5000)
    fpres, ftotal = rep(fardata, 'FarfromReg')
    cpres, ctotal = rep(closedata, 'CloseToReg') 
    return fardata, closedata
    

def meta(file = 'GivenFiles/CMAIKI_Metadata.csv'):
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

def addMD():
    md = meta()
    data = imp()
    df = pd.concat([data, md], sort = False, join = 'inner', axis = 0)
    return df

def normMD():
    md = meta()
    data = imp()
    data = data.div(data.sum(axis = 0))
    df = pd.concat([data, md], sort = False, join = 'inner', axis = 0)
    return df

def qualanalysis(key, df = addMD()):
    cats = [];
    for cat in df.loc[key]:
        if cat in cats:
            continue
        else:
            cats.append(cat)
    vals = []
    l = int(len(cats))
    n = math.sqrt(l)
    n = math.ceil(n)
    for cat in cats:
        val = df.loc[key] == cat
        vals.append(val)
    df = df.transpose()
    x = "Presence"
    y = "Total"
    markers = ['b.', 'gv', 'r1', 'cs', 'mp', 'y*', 'kx', 'bv', 'g1', 'rs', 'cp', 'm*', 'yx', 'k.',
               'b1', 'gs', 'rp', 'c*', 'mx', 'y.', 'kv', 'bs', 'gp', 'r*', 'cx', 'm.', 'yv', 'k1', 
               'bp', 'g*', 'rx', 'c.', 'mv', 'y1', 'ks', 'b*', 'gx', 'r.', 'cv', 'm1', 'ys', 'kp']
    fig = plt.figure(figsize = (4*n,3*n))
    for i in range(l):
        val = vals[i]
        tempdf = df[val]
        tempdf = tempdf.transpose()
        per = round((tempdf.shape[1]/df.shape[0])*100)
        title = str(per) + "% of samples are " +str(cats[i]) 
        p, t = summarize(tempdf.iloc[0:-4, ])
        fn, residual , rank, singularValues, condition = np.polyfit(p, t, 1, full = True)
        f =  np.poly1d(fn)
        fig.add_subplot(n, n, i+1)
        plt.plot(p, t, markers[i])
        plt.plot(p, f(p), 'k-', label = 'unfiltered m = ' + str(round(fn[0], 2)))
        lim = math.sqrt(residual[0])/ math.sqrt(tempdf.shape[1])
        plt.plot(p, [lim]*len(p), 'g-', label = 'filtering cutoff')
        tempdf = filt(tempdf.iloc[0:-4, ], t, lim)
        fp, ft = summarize(tempdf)
        fn = np.polyfit(fp, ft, 1)
        f =  np.poly1d(fn)
        plt.plot(p, f(p), 'r-', label = 'filtered m = '+str(round(fn[0], 2)))
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2))
        plt.title(title)
        plt.xlabel(x)
        plt.ylabel(y)
    fig.tight_layout()
    fig.savefig(key+"Plots.jpg")
    fig = plt.figure(figsize = (8, 8))
    for i in range(l):
        val = vals[i]
        tempdf = df[val]
        tempdf = tempdf.transpose()
        title = key 
        p, t = summarize(tempdf.iloc[0:-4, ])
        plt.plot(p, t, markers[i], label = cats[i])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(title)
    plt.xlabel(x)
    plt.ylabel(y)
    fig.tight_layout()
    fig.savefig(key+"Plot.jpg")
    
def qualAnal():
    #df = addMD()
    cols = [ 'SampleType', 'Habitat', 'Host', 'Trophic']
    df = normMD()
    for key in cols:
        qualanalysis(key, df)
        
def connect(DM, e): # adjacency matrix
    C = [];
    for x in DM:
        if x >= e: C.append(1)
        else: C.append(0)
    return C

def adjacency(DM, e): # adjacency matrix
    C = [];
    for x in DM:
        if x >= e: C.append(x) # keep connections that are similar in exp(-d^2) form 1 is similar, 0 is far
        else: C.append(0)
    return C

def degree(C, l):
    tri = CtoSymM(C,l)
    d = np.sum(tri, axis = 1)
    return d

def CtoSymM(C, l):
    tri = np.zeros((l,l))
    tri[np.triu_indices(l, 1)] = C
    tri = tri + tri.transpose()
    return tri


def bcdist(u,v):
    n = len(u)
    c = 0;
    for i in range(n):
        if u[i]*v[i] != 0:
            c += min(u[i], v[i])
    return c

def BrayCurtis(df):
    n = df.shape[1]
    dm = []
    for i in range(n):
        for j in range(i+1, n):
            u = df.iloc[:,i]
            v = df.iloc[:,j]
            c = bcdist(u,v)
            dm.append(c)
    return np.array(dm)

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
    


#def main():
#
#
#    
#            
#if __name__ == "__main__":
#    main()
#    
 

