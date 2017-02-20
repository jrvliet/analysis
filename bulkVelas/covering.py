#!/usr/bin/python

'''
Plots covering fraction over time as a heatmap
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

# Cf parameters
CL = 0.8413
tol = 13-4
pmin = 0.0
pmax = 1.0
iterations = 1000

rootloc = '/mnt/cluster/abs/cgm/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i{2:d}/{3:s}/'
filename = '{0:s}.{1:s}.a{2:.3f}.ALL.sysabs'

galNums = range(20,30)

ions = 'HI MgII CIV OVI'.split()
ewcut = 0.1

loD,hiD = 0,1.5
numDbins = 15
Dbins = np.linspace(loD,hiD,numDbins+1)
DbinLabels = ['{0:.1f}'.format(i) for i in Dbins[1:]]

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
        galID,expn,redshift,mvir,rvir,inc = galaxyProps(loc)
        
        for ion in ions:

            loc = rootloc+subloc.format(galNum,a,inc,ion)
            sysabs = loc+filename.format(galID,ion,a)
            impact,ew = np.loadtxt(sysabs,skiprows=1,usecols=(1,5),unpack=True)
        
            linesfile = loc+'lines.info'
            lines = np.loadtxt(linesfile,skiprows=2)

            # The impact parameters need to be rounded becuase the impact parameters
            # in the ALL.sysabs file are rounded to 1 decimal point
            linesImp = np.round(lines[:,1],1)

            for i in range(numDbins):
                loD = np.round(Dbins[i]*rvir,1)
                hiD = np.round(Dbins[i+1]*rvir,1)
                ews = ew[(impact>=loD) & (impact<hiD)]

                numHits = (ews>ewcut).sum()
                numLOS = ((linesImp>=loD) & (linesImp<hiD)).sum()
                fraction = float(numHits)/float(numLOS)
                results[aLabel,ion].iloc[i] = fraction

    s = 'vela2b-{0:d}_covering.h5'.format(galNum)
    results.to_hdf(s,'data',mode='w')
       



