

from __future__ import print_function
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None

loc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'

expns = np.arange(0.200,0.550,0.01)

lowestN,highestN = -5.5,-2.5
numNbins = 3
nBins = np.linspace(lowestN,highestN,numNbins+1)

lowestT,highestT = 3.5,6.5
numTbins = 3
tBins = np.linspace(lowestT,highestT,numTbins+1)

denseLabels = 'diffuse mid dense'.split()
tempLabels = 'cool warm hot'.split()
fields = 'radVelRatio fracInfall'.split()
expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]

header = [tempLabels,denseLabels,fields]
header = pd.MultiIndex.from_product(header)
results = np.zeros((len(expns),len(header)))
results = pd.DataFrame(results,columns=header,index=expnLabels)

for a,aLabel in zip(expns,expnLabels):

    print(aLabel)
    fname = loc+filename.format(a)
    df = pd.read_hdf(fname,'data')

    df['r'] = np.sqrt(df['xRot']**2 + df['yRot']**2 + df['zRot']**2)
    df['vr'] = (df['vxRot']*df['xRot'] + df['vxRot']*df['xRot'] + 
                df['vxRot']*df['xRot']) / df['r']

    df['speed'] = np.sqrt(df['vxRot']**2 + df['vyRot']**2 + df['vzRot']**2)

    spaceInds = df['theta']<80

    for i,tempLabel in enumerate(tempLabels):

        loT,hiT = tBins[i],tBins[i+1]
        tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
        
        for j,denseLabel in enumerate(denseLabels):

            loN,hiN = nBins[j],nBins[j+1]
            denseInds = (df['density']>10**loN) & (df['density']<10**hiN)

            fil = df[spaceInds & tempInds & denseInds]
            
            rVelRatio = fil['vr']/fil['speed']
            fracInfall = float((fil['vr']<0).sum()) / len(fil)

            values = [rVelRatio.mean(), fracInfall]
            for field,val in zip(fields,values):
                results[tempLabel,denseLabel,field][aLabel] = val
            
hdfName = 'fractionInfall.h5'
results.to_hdf(hdfName,'data')
    

