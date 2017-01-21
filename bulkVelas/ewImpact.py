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
subloc = 'a{0:.3f}/i{1:d}/{2:s}/'
filename = '{0:s}.{1:s}.a{2:.3f}.ALL.sysabs'
expns = np.arange(0.200,0.540,0.01)
expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]

ions = 'HI MgII CIV OVI'.split()

loD,hiD = 0,1.5
numDbins = 15
Dbins = np.linspace(loD,hiD,numDbins+1)
DbinLabels = ['{0:.1f}'.format(i) for i in Dbins[1:]]

header = [expnLabels,ions]
header = pd.MultiIndex.from_product(header)
results = np.zeros((numDbins,len(header)))
results = pd.DataFrame(results,columns=header,index=DbinLabels)

for a,aLabel in zip(expns,expnLabels):

    print(a)
    loc = rootloc+'a{0:.3f}/'.format(a)
    galID,expn,redshift,mvir,rvir,inc = galaxyProps(loc)
    
    for ion in ions:

        print('\t',ion)
        sysabs = rootloc+subloc.format(a,inc,ion)+filename.format(galID,ion,a)
        impact,ew = np.loadtxt(sysabs,skiprows=1,usecols=(1,5),unpack=True)

        
        for i in range(numDbins):
            loD = Dbins[i]*rvir
            hiD = Dbins[i+1]*rvir
            ews = ew[(impact>=loD) & (impact<hiD)]
            results[aLabel,ion].iloc[i] = np.mean(ews)

s = 'ewImpact.h5'
results.to_hdf(s,'data',mode='w')
        


        
    



