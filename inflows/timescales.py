
'''
Determines mean/std of gas timescales in the filament
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


loc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.{1:s}.h5'

expns = np.arange(0.200,0.550,0.01)
aLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]

ions = 'HI MgII CIV OVI'.split()
fields = 't_ph t_rec t_coll t_cool'.split()

lowestT,highestT = 3.5,6.5
lowestN,highestN = -5.5,-2.5
numTbins,numNbins = 3,3
tBins = np.linspace(lowestT,highestT,numTbins+1)
nBins = np.linspace(lowestN,highestN,numNbins+1)
tLabels = 'cool warm hot'.split()
nLabels = 'diffuse mid dense'.split()

header = [tLabels,nLabels,fields]
header = pd.MultiIndex.from_product(header)

for ion in ions:
    print(ion)

    times = pd.DataFrame(columns=header,index=aLabels)

    for a,aLabel in zip(expns,aLabels):

        print('\t',a)

        fname = loc+filename.format(a,ion)
        df = pd.read_hdf(fname,'data')
        rotname = loc+filename.format(a,'rot')
        rotdf = pd.read_hdf(rotname,'data')

        spaceInds = rotdf['theta']<80

        for i,tLabel in enumerate(tLabels):

            loT,hiT = tBins[i],tBins[i+1]
            tempInds = (df['temperature']<10**hiT) & (df['temperature']>10**loT)
            
            for j,nLabel in enumerate(nLabels):

                loN,hiN = nBins[j],nBins[j+1]
                denseInds = (df['nH']<10**hiN) & (df['nH']>10**loN)

                fil = df[spaceInds & tempInds & denseInds]
                for field in fields:
                    val = np.log10((10**fil[field]).mean())
                    times[tLabel,nLabel,field].ix[aLabel] = val
    
               
    s = 'timescales_{0:s}.h5'.format(ion)
    times.to_hdf(s,'data',mode='w')




