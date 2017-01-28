#!/usr/bin/python

'''
Plots binned EW vs. D over time as a heatmap
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None

def galaxyProps(location):

    fname = location+'galaxy.props'
    with open(fname) as f:
        galID = f.readline().split()[1]
        expn = f.readline().split()[1]
        redshift = f.readline().split()[1]
        mvir = f.readline().split()[1]
        rvir = float(f.readline().split()[1])
        inc = int(float(f.readline().strip().split()[1]))

    return galID,expn,redshift,mvir,rvir,inc

rootloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
rootloc = '/home/jacob/research/velas/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i{2:d}/{3:s}/'
filename = '{0:s}.{1:s}.a{2:.3f}.i{3:d}.ALL.sysabs.h5'

ions = 'HI MgII CIV OVI'.split()
ions = 'HI MgII CIV'.split()

loD,hiD = 0,1.5
numDbins = 15
Dbins = np.linspace(loD,hiD,numDbins+1)
DbinLabels = ['{0:.1f}'.format(i) for i in Dbins[1:]]

finalExpn = [0.530, 0.490, 0.490, 0.510, 0.460, 
             0.500, 0.500, 0.500, 0.510, 0.500]

galNums = range(20,30)
for galNum,finala in zip(galNums,finalExpn):

    print(galNum)
    expns = np.arange(0.200,finala,0.01)
    expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]

    header = [expnLabels,ions]
    header = pd.MultiIndex.from_product(header)
    results = np.zeros((numDbins,len(header)))
    results = pd.DataFrame(results,columns=header,index=DbinLabels)

    for a,aLabel in zip(expns,expnLabels):

        print(a)
        loc = rootloc+'vela{0:d}/a{1:.3f}/i90/'.format(galNum,a)
        galID,expn,redshift,mvir,rvir,inc = galaxyProps(loc)
        
        print(type(inc))
        for ion in ions:

            print('\t',ion)
            sysabs = rootloc+subloc.format(galNum,a,inc,
                        ion)+filename.format(galID,ion,a,inc)
            df = pd.read_hdf(sysabs,'data')

            for i in range(numDbins):
                loD = Dbins[i]*rvir
                hiD = Dbins[i+1]*rvir
                ews = df['EW_r'][(df['D']>=loD) & (df['D']<hiD)]
                results[aLabel,ion].iloc[i] = ews.mean()

    s = 'vela2b-{0:d}_ewImpact.h5'.format(galNum)
    results.to_hdf(s,'data',mode='w')

    break

