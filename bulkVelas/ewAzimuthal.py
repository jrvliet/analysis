#!/usr/bin/python

'''
Plots binned EW vs. azimuthal angle over time as a heatmap
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None
import sys

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

def azimuthal(phi):

    if phi<90:
        alpha = phi
    elif phi<180:
        alpha = 180-phi
    elif phi<270:
        alpha = phi-180
    else:
        alpha = 360-phi

    return alpha

rootloc = '/mnt/cluster/abs/cgm/vela2b/'
rootloc = '/home/jacob/research/velas/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i{2:d}/{3:s}/'
filename = '{0:s}.{1:s}.a{2:.3f}.i{3:d}.ALL.sysabs.h5'

galNums = range(21,30)

ions = 'HI MgII CIV OVI'.split()

loAz,hiAz = 0,90
numAzbins = 10
azBins = np.linspace(loAz,hiAz,numAzbins+1)
azBinLabels = ['{0:.1f}'.format(i) for i in azBins[1:]]

finalExpn = [0.530, 0.490, 0.490, 0.510, 0.460, 
             0.500, 0.500, 0.500, 0.510, 0.500]

for galNum,finala in zip(galNums,finalExpn):

    print(galNum)
    expns = np.arange(0.200,finala,0.01)
    expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]

    header = [expnLabels,ions]
    header = pd.MultiIndex.from_product(header)
    results = np.zeros((numAzbins,len(header)))
    results = pd.DataFrame(results,columns=header,index=azBinLabels)

    for a,aLabel in zip(expns,expnLabels):

        print(a)
        try:
            loc = rootloc+'vela{0:d}/a{1:.3f}/i90/'.format(galNum,a)
            galID,expn,redshift,mvir,rvir,inc = galaxyProps(loc)
            
            for ion in ions:

                print('\t',ion)
                sysabs = rootloc+subloc.format(galNum,a,inc,
                            ion)+filename.format(galID,ion,a,inc)
                df = pd.read_hdf(sysabs,'data')
                df['azimuthal'] = df['phi'].apply(azimuthal)

                for i in range(numAzbins):
                    loAz = azBins[i]
                    hiAz = azBins[i+1]
                    ews = df['EW_r'][(df['azimuthal']>=loAz) & (df['azimuthal']<hiAz)]
                    results[aLabel,ion].iloc[i] = ews.mean()
        except IOError:
            continue

    s = 'vela2b-{0:d}_ewAzimuthal.h5'.format(galNum)
    results.to_hdf(s,'data',mode='w')


