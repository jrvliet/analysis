
from __future__ import print_function
#import matplotlib as mp
#mp.use('Agg')
#import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#import scipy.stats as st
#import matplotlib.ticker as ticker
pd.options.mode.chained_assignment = None
#import scipy.optimize as opt

kpc2km = 3.086e16 
pc2cm = 3.086e18 
mH = 1.6737e-24 
mSun = 1.989e33 
s2yr = 3.154e7 

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

# Expansion Parameter steps
expns = np.arange(0.200,0.550,0.01)
expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]

# Density bins
lowestN,highestN = -5.5,-2.5
numNbins = 3
denseBins = np.linspace(lowestN,highestN,numNbins+1)
denseLabels = []
for i in range(numNbins):
    s = 'nLim=({0:.1f},{1:.1f})'.format(denseBins[i],denseBins[i+1])
    denseLabels.append(s)
denseLabels = 'rare mid dense'.split()

# Temperatre bins
lowestT,highestT = 3.5,6.5
numTbins = 3
tempBins = np.linspace(lowestT,highestT,numTbins+1)
tempLabels = 'cool warm hot'.split()

# Set up dataframe
fields = 'm90 rho90Rvir rho90kpc'.split()
header = [tempLabels,denseLabels,fields]
header = pd.MultiIndex.from_product(header)
results = np.zeros((len(expns),len(header)))
results = pd.DataFrame(results,columns=header,index=expnLabels)

for i,a in enumerate(expns):

    print(a)
    expnLabel = expnLabels[i]

    rname = dataloc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = dataloc + filename.format(a)
    df = pd.read_hdf(fname,'data')

    # Select out full filament
    fullTempInds = ((df['temperature']>10**lowestT) &
                    (df['temperature']<10**highestT))
    spaceInds = (df['theta']<80) & (df['r']>0.5*rvir)
    fullDenseInds = (df['density']>10**lowestN) & (df['density']<10**highestN)
    fullFilament = df[fullTempInds & spaceInds & fullDenseInds]

    fullFilament['mass'] = fullFilament['density']*mH*(fullFilament['cell_size']*pc2cm)**3 / mSun
    print('Fraction of box in filament = {0:f}'.format( 
                    len(fullFilament)/float(len(df))))
    totalMass = fullFilament['mass'].sum()

    # Loop over temperature bins
    for j in range(numTbins):
        tempLabel = tempLabels[j]
        print('\t{0:s}'.format(tempLabel))
        loT = tempBins[j]
        hiT = tempBins[j+1]
        tempInds = ((fullFilament['temperature']>10**loT) &
                    (fullFilament['temperature']<10**hiT))

        # Loop over density bins
        for k in range(numNbins):
        
            denseLabel = denseLabels[k]
            print('\t\t{0:s}'.format(denseLabel))
            loN = denseBins[k]
            hiN = denseBins[k+1]
            denseInds = ((fullFilament['density']>10**loN) & 
                            (fullFilament['density']<10**hiN))

            fil = fullFilament[tempInds & denseInds]
            fil['rho'] = np.sqrt(fil['xRot']**2 + fil['yRot']**2)/rvir

            rhoMin,rhoMax = 0,3
            rhoBins = np.linspace(rhoMin,rhoMax,5000)

            # Calculate mass contained within each rho bin
            mIn = []
            for k in range(len(rhoBins)):
                inside = fil['rho']<=rhoBins[k]
                massInside = fil['mass'][inside].sum()
                mIn.append(massInside)

            # Find the value of rho which contains 90% of mass
            massFraction = 0.90
            m90 = mIn[-1]*massFraction
            percentile = np.digitize(np.array(m90),mIn)
            if percentile==len(rhoBins):
                containingRadius = rhoBins[-1]
            else:
                containingRadius = rhoBins[percentile]

            values = [m90,containingRadius,containingRadius*rvir]
            for field,val in zip(fields,values):
                results[tempLabel,denseLabel,field][expnLabel] = val

hdfName = 'filamentMassSize_simple.h5'
results.to_hdf(hdfName,'data')

