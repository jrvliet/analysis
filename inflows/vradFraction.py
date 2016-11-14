
'''
Plots the fraction of the velocity in the filament that is in the radial
direction for various density cuts
'''


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

pd.options.mode.chained_assignment = None


dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'

expns = np.arange(0.200,0.500,0.01)

loT,hiT = 3.5,4.5
lowestN,highestN = -5.5,-2.5
numDenseBins = 6
denseEdges = np.linspace(lowestN,highestN,numDenseBins+1)

labels = []
for i in range(numDenseBins):
    labels.append('{0:.1f}<n<{1:.1f}'.format(denseEdges[i],denseEdges[i+1]))


fraction = np.zeros((len(expns),numDenseBins))
spread = np.zeros((len(expns),numDenseBins))

for i,a in enumerate(expns):

    print(a)

    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname,'data')

    df['speed'] = np.sqrt(df['vxRot']**2 + df['vyRot']**2 + df['vzRot']**2)
    df['vr'] = (df['vxRot']*df['xRot'] + df['vyRot']*df['yRot'] +
                df['vzRot']*df['zRot']) / df['r']

    tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    spacInds = (df['theta']<80)
    
    for j in range(numDenseBins):

        print('\t{0:s}'.format(labels[j]))
        loN = denseEdges[j]
        hiN = denseEdges[j+1]

        densInds = (df['density']<10**hiN) & (df['density']>10**loN)
        
        cloud = df[densInds & tempInds & spacInds]

        cloud['vrFrac'] = cloud['vr'] / cloud['speed']
        
        fraction[i,j] = cloud['vrFrac'].mean()
        spread[i,j] = cloud['vrFrac'].std()

   
dfFrac = pd.DataFrame(fraction,columns=labels,index=expns)
dfSpread = pd.DataFrame(spread,columns=labels,index=expns)

dfFrac.to_hdf('vradFraction_mean.h5','data',mode='w')
dfSpread.to_hdf('vradFraction_std.h5','data',mode='w')



        

