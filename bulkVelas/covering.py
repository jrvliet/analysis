#!/usr/bin/python

'''
Plots covering fraction over time as a heatmap
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None

def galaxyProps(location):

    fname = location+'/i90/galaxy.props'
    with open(fname) as f:
        galID = f.readline().split()[1]
        expn = f.readline().split()[1]
        redshift = f.readline().split()[1]
        mvir = f.readline().split()[1]
        rvir = float(f.readline().split()[1])
        inc = int(float(f.readline().strip().split()[1]))

    return galID,expn,redshift,mvir,rvir,inc

# Cf parameters
CL = 0.8413
tol = 13-4
pmin = 0.0
pmax = 1.0
iterations = 1000

rootloc = '/mnt/cluster/abs/cgm/vela2b/'
rootloc = '/home/jacob/research/velas/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i{2:d}/{3:s}/'
filename = '{0:s}.{1:s}.a{2:.3f}.i90.ALL.sysabs.h5'

galNums = range(21,30)

ions = 'HI MgII CIV OVI'.split()
ewcut = 0.1

loD,hiD = 0,1.5
numDbins = 15
Dbins = np.linspace(loD,hiD,numDbins+1)
DbinLabels = ['{0:.1f}'.format(i) for i in Dbins[1:]]

finalExpn = [0.550]*len(galNums)

for galNum,finala in zip(galNums,finalExpn):

    print(galNum)

    expns = np.arange(0.200,finala,0.01)
    expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]
    header = [expnLabels,ions]
    header = pd.MultiIndex.from_product(header)
    results = np.zeros((numDbins,len(header)))
    results = pd.DataFrame(results,columns=header,index=DbinLabels)

    for a,aLabel in zip(expns,expnLabels):

        loc = rootloc+'vela{0:d}/a{1:.3f}/'.format(galNum,a)

        try:
            galID,expn,redshift,mvir,rvir,inc = galaxyProps(loc)
            
            for ion in ions:

                loc = rootloc+subloc.format(galNum,a,inc,ion)
                sysabs = loc+filename.format(galID,ion,a)
                df = pd.read_hdf(sysabs,'data')
            
                for i in range(numDbins):
                    loD = np.round(Dbins[i]*rvir,1)
                    hiD = np.round(Dbins[i+1]*rvir,1)

                    index = (df['D']>=loD) & (df['D']<hiD)
                    d = df[index]

                    numHits = (d['EW_r']>ewcut).sum()
                    numLOS = index.sum()

                    fraction = float(numHits)/float(numLOS)
                    results[aLabel,ion].iloc[i] = fraction

        except IOError:
                continue

    s = 'vela2b-{0:d}_covering.h5'.format(galNum)
    results.to_hdf(s,'data',mode='w')
       



