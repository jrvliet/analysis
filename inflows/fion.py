
'''
Plots the evolution of fion in each filament phase
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None

u2g = 1.661e-24
g2M = 1.989e33
pc2cm = 3.086e18
boltz = 1.3807e-16
cm2km = 1e5

loc = '/home/jacob/research/velas/vela2b/vela27/'
ionfilename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.{1:s}.h5'
boxfilename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
ions = 'HI MgII CIV OVI'.split()
expns = np.arange(0.200,0.500,0.01)
expns = [0.400]

lowestT,highestT = 3.5,6.5
lowestN,highestN = -2.5,-5.5
numTbins,numNbins = 3,3
tBins = np.linspace(lowestT,highestT,numTbins+1)
nBins = np.linspace(lowestN,highestN,numNbins+1)

tempHeader = 'cool warm hot'.split()
denseHeader = 'diffuse mid dense'.split()
fields = 'mean std'.split()
header = [tempHeader,denseHeader,fields]
header = pd.MultiIndex.from_product(header)
expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]
results = pd.DataFrame(columns=header,index=expnLabels)

for i,a in enumerate(expns):
    print(a)

    boxf = loc+boxfilename.format(a)
    box = pd.read_hdf(boxf,'data')

    spaceInds = box['theta']<80

    print(spaceInds.sum())
    for ion in ions:

        print(ion)

        fname = loc+ionfilename.format(a,ion)
        df = pd.read_hdf(fname,'data')

        for j in range(numTbins):
            
            tempLabel = tempHeader[j]
            loT,hiT = tBins[j],tBins[j+1]
            tempBins = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
            print(tempBins.sum())
            
            for k in range(numNbins):

                denseLabel = denseHeader[k]
                loN,hiN = nBins[k],nBins[k+1]
                denseBins1 = (df['nH']<10**loN) 
                denseBins2 = (df['nH']>10**hiN)
                denseBins = denseBins1 & denseBins2
                print(denseBins1.sum())
                print(denseBins2.sum())
                print(loN,hiN)
                print(df['nH'].describe())
                print(denseBins.sum())

                fil = df[spaceInds & tempBins & denseBins]
                print(len(fil))

                fil['mass'] = (fil['nH']*u2g*pc2cm**3/g2M)*fil['cell_size']**3
    
                results[tempLabel,denseLabel,fields[0]].ix[i] = np.average(fil['fIon'],
                                                                    weights=fil['mass'])
                results[tempLabel,denseLabel,fields[1]].ix[i] = fil['fIon'].std()

s = 'fIon.h5'
results.to_hdf(s,'data',mode='w')
                
                
                 
    


