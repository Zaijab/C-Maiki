#!/usr/bin/env python

# Import necessary Modules
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline
from datetime import datetime
from itertools import cycle
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Import Data from Excell as Pandas Dataframe
def createdf(path):
    tempSheet = pd.read_excel(path,)
    dates = pd.date_range(start='2015-08-05',end='2016-07-10',freq='30T')
    tempsdict = {"10m": tempSheet['Unnamed: 2'][13:16334].to_numpy(),'30m': tempSheet['Unnamed: 4'][157:-2].to_numpy(), '40m': tempSheet['Unnamed: 6'][13:16334].to_numpy(),'60m': tempSheet['Unnamed: 8'][13:16334].to_numpy(),'70m': tempSheet['Unnamed: 10'][13:16334].to_numpy(),'80m': tempSheet['Unnamed: 12'][13:16334].to_numpy(),'90m': tempSheet['Unnamed: 14'][13:16334].to_numpy(),'120m': tempSheet['Unnamed: 16'][13:16334].to_numpy(),'130m': tempSheet['Unnamed: 18'][13:16334].to_numpy()}
    return pd.DataFrame(data=tempsdict, index=dates)

# Define Distances
def c1_norm(function1,function2):
    # Take the difference of the functions
    functiondiff = np.array(function1) - np.array(function2)
    interpolant = CubicSpline([i for i in range(len(functiondiff))],functiondiff,bc_type='natural').derivative().__call__([i for i in range(len(functiondiff))])
    return np.sqrt(sum(np.square(functiondiff)+np.square(interpolant)))

def old_c1(function1,function2):
    return np.sqrt(c1_no_derivative(function1,function2)+c1_derivative(function1,function2))

# Interesting Parts of Distances
def c1_no_derivative(function1,function2):
    return (sum(np.square(np.asarray(function1) - np.asarray(function2))))

def c1_derivative(function1,function2):
    localchange = []
    functiondiff = cycle(np.asarray(function1) - np.asarray(function2))
    for index in range(len(function1)):
        locallist = []
        for localindex in range(3):
            locallist.append(next((functiondiff)))


        for i in range(len(function1)-2):
            next((functiondiff))

        localchange.append(np.polyfit([1,2,3],locallist,1)[0])
        return sum(np.square(localchange))

def c1_spline_derivative(function):
    return np.sqrt(sum(np.square(CubicSpline([i for i in range(len(function))],function,bc_type='natural').derivative().__call__([i for i in range(len(function))]))))

# Plotting Helper Functions
def plotter(matrix,header=None,title=None,clim=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.imshow(matrix, interpolation='nearest')
    fig.colorbar(cax)
    #if clim != None:
    #cax.set_clim(0,100)
    if header != None:
        ax.set_xticklabels(['']+header)
        ax.set_yticklabels(['']+header)
    if title != None:
        ax.set_title(title)
    ax.xaxis.set_major_locator(ticker.FixedLocator([0,1,2,3,4,5,6,7,8]))
    ax.yaxis.set_major_locator(ticker.FixedLocator([0,1,2,3,4,5,6,7,8]))
    #ax.xaxis.set_major_locator(ticker.FixedLocator([0,31,61,92,122,153,184,213,244,274,305,335]))
    #ax.yaxis.set_major_locator(ticker.FixedLocator([0,31,61,92,122,153,184,213,244,274,305,335]))
    ax.set_yticklabels(["10m","30m","40m","60m","70m","80m","90m","120m","130m"])
    ax.set_xticklabels(["10m","30m","40m","60m","70m","80m","90m","120m","130m"])
    #ax.set_xticklabels(["August","September","October","November","December","January","February","March","April","May","June","July"])
    #ax.set_yticklabels(["August","September","October","November","December","January","February","March","April","May","June","July"])

    plt.xticks(rotation='vertical')
    plt.show()

# Calculate Distance Matrices
def distance_matrix(vectors):
    distmatrix = np.zeros((len(vectors),len(vectors)))
    for i in range(len(vectors)):
        for j in range(len(vectors)):
            distmatrix[i][j] = c1_norm((vectors[i]),(vectors[j]))
    return distmatrix

