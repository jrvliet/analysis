
'''
Plots the distribution of pressure in gas when binned by temperature
(cool, warm, hot) and metallicity
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def mfToZ(mf):

    X_sun = 0.70683
    Y_sun = 0.27431
    Z_sun = 0.0188 
    r = 0.3347
    #X_cell = (1 - mf) / (1 + r)
    #Z = 
    return (mf / ((1 - mf) / (1 + r))) / (Z_sun / X_sun)


dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expns = np.arange(0.200,0.550,0.01)

boltz = 1.380658e-16

loN,hiN = -5.5,-2.5
loT,hiT = 3.5,4.5

header = 'a min 25 50 75 max mean'.split()
header = ['a loMean loStd '
            'loMidMean loMidStd ' 
            'hiMidMean hiMidStd '
            'hiMean hiStd']
header = header[0].split()
statsArr = np.zeros((len(expns),len(header)))
stats = pd.DataFrame(statsArr,columns=header)

zBins = [-12.5,-2.62,-2.38,-2.12,-1.9]
zLabels = 'lo loMid hiMid hi'.split()

for i,a in enumerate(expns):
    
    print(a)
    
    rname = dataloc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = dataloc+filename.format(a)

    df = pd.read_hdf(fname,'data')

    df['r'] = np.sqrt(df['xRot']**2 + df['yRot']**2 + df['zRot']**2 )/rvir

    tempInd = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    densInd = (df['density']>10**loN) & (df['density']<10**hiN)
    spacInd = (df['theta']<80)
    distInd = (df['r']>0.5)
    
    df = df[tempInd & densInd & spacInd & distInd]

    df['metal'] = np.log10(mfToZ(df['SNII']+df['SNIa']))
    df['pressure'] = df['density']*boltz*df['temperature']
    #df['metal'] = np.log10(mfToZ(df['SNII']))

    fig,ax = plt.subplots(1,1,figsize=(5,5))
    for j in range(len(zBins)-1):
    
        loZ = zBins[j]
        hiZ = zBins[j+1]
        zBinLabel = zLabels[j]

        zInd = (df['metal']>loZ) & (df['metal']<hiZ)
        cl = df[zInd]

        print('\t{0:d}'.format(zInd.sum()))

        stats[zBinLabel+'Mean'].ix[i] = cl['pressure'].mean()
        stats[zBinLabel+'Std'].ix[i] = cl['pressure'].std()

        ax.hist(np.log10(cl['pressure']),bins=25,log=True,histtype='step',label=zBinLabel)
        
    ax.legend(loc='best')
    
    ax.set_title(a)
    fig.savefig('pressure_dist_a{0:.3f}.png'.format(a))
    plt.close(fig)

s = 'pressureStats_noCGM.h5'
stats.to_hdf(s,'data',mode='w')



