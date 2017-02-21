#!/usr/bin/python

'''
Plots covering fraction vs azimuthal over time as a heatmap
'''

from __future__ import print_function
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



# Cf parameters
CL = 0.8413
tol = 13-4
pmin = 0.0
pmax = 1.0
iterations = 1000

rootloc = '/mnt/cluster/abs/cgm/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i{2:d}/{3:s}/'
filename = '{0:s}.{1:s}.a{2:.3f}.i90.ALL.sysabs.h5'

galNums = range(20,30)
galNums = range(21,29)

finalExpn = [0.510, 0.490, 0.490, 0.510, 0.460,
             0.500, 0.500, 0.500, 0.510, 0.500]



ions = 'HI MgII CIV OVI'.split()
ewcut = 0.1

loAz,hiAz = 0,90
numAzbins = 10
Azbins = np.linspace(loAz,hiAz,numAzbins+1)
AzbinLabels = ['{0:d}'.format(int(i)) for i in Azbins[1:]]

print(Azbins)
print(AzbinLabels)

for galNum,finala in zip(galNums,finalExpn):

    print(galNum)

    expns = np.arange(0.200,finala,0.01)
    expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]
    fields = 'full lowD highD'.split()
    header = [expnLabels,ions,fields]
    header = pd.MultiIndex.from_product(header)
    results = np.zeros((numAzbins,len(header)))
    results = pd.DataFrame(results,columns=header,index=AzbinLabels)

    for a,aLabel in zip(expns,expnLabels):

        print('\t',aLabel)
        loc = rootloc+'vela{0:d}/a{1:.3f}/'.format(galNum,a)
        galID,expn,redshift,mvir,rvir,inc = galaxyProps(loc)
        
        for ion in ions:

            loc = rootloc+subloc.format(galNum,a,inc,ion)
            sysabs = loc+filename.format(galID,ion,a)
            df = pd.read_hdf(sysabs,'data')
            df['azi'] = df['phi'].apply(azimuthal)
        
            linesfile = loc+'lines.info'
            lines = np.loadtxt(linesfile,skiprows=2)

            # The impact parameters need to be rounded becuase the impact parameters
            # in the ALL.sysabs file are rounded to 1 decimal point
            linesImp = np.round(lines[:,1],1)

            for i in range(numAzbins):
                loAz,hiAz = Azbins[i],Azbins[i+1]

                azIndex = (df['azi']>=loAz) & (df['azi']<hiAz)
                loDIndex = (df['D']<0.3*rvir) & azIndex
                hiDIndex = (df['D']>=0.3*rvir) & azIndex

                indicies = [azIndex,loDIndex,hiDIndex]
                for field,index in zip(fields,indicies):
                    d = df[index]
                    numHits = (d['EW_r']>ewcut).sum()
                    numLOS = index.sum()
                    if numLOS>0:
                        fraction = float(numHits)/float(numLOS)
                    else:
                        fraction = 0
                    results[aLabel,ion,field].iloc[i] = fraction
                    

       
    s = 'vela2b-{0:d}_covering.h5'.format(galNum)
    results.to_hdf(s,'data',mode='w')
    break
       



