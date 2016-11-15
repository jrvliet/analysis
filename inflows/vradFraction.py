
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
import itertools as it

pd.options.mode.chained_assignment = None

# Filenames
dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

# Snapshots
expns = np.arange(0.200,0.500,0.01)

# Temperature bins
numTempBins = 3
lowestT, highestT = 3.5, 6.5
tempEdges = np.linspace(lowestT,highestT,numTempBins+1)
tempLabels = 'cool warm hot'.split()

# Density bins
loT,hiT = 3.5,4.5
lowestN,highestN = -5.5,-2.5
numDenseBins = 6
denseEdges = np.linspace(lowestN,highestN,numDenseBins+1)

denseLabels = []
for i in range(numDenseBins):
    denseLabels.append('{0:.1f}<n<{1:.1f}'.format(denseEdges[i],denseEdges[i+1]))

# Initialize zero arrays for output
vrFraction = np.zeros((numTempBins,len(expns),numDenseBins))
vrStd = np.zeros((numTempBins,len(expns),numDenseBins))
rMean = np.zeros((numTempBins,len(expns),numDenseBins))
rStd = np.zeros((numTempBins,len(expns),numDenseBins))

# Loop over snapshots
for i,a in enumerate(expns):

    print(a)

    rname = dataloc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname,'data')

    df['speed'] = np.sqrt(df['vxRot']**2 + df['vyRot']**2 + df['vzRot']**2)
    df['vr'] = (df['vxRot']*df['xRot'] + df['vyRot']*df['yRot'] +
                df['vzRot']*df['zRot']) / df['r']

    spacInds = (df['theta']<80)
    # Loop over temperature
    for j in range(numTempBins):
        print('\t{0:s}'.format(tempLabels[j]))
        loT = tempEdges[j]
        hiT = tempEdges[j+1]
        tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    
        for k in range(numDenseBins):

            print('\t\t{0:s}'.format(denseLabels[k]))
            loN = denseEdges[k]
            hiN = denseEdges[k+1]

            densInds = (df['density']<10**hiN) & (df['density']>10**loN)
            
            cloud = df[densInds & tempInds & spacInds]

            cloud['vrFrac'] = cloud['vr'] / cloud['speed']
            
            vrFraction[j,i,k] = cloud['vrFrac'].mean()
            vrStd[j,i,k] = cloud['vrFrac'].std()
            rMean[j,i,k] = cloud['r'].mean()/rvir
            rStd[j,i,k] = cloud['r'].std()/rvir

   
# Ouput files
outputFields = 'vradFrac r'.split()
outputStats = 'mean std'.split()
stores = []
arrs = [vrFraction, vrStd, rMean, rStd]
for i,combos in enumerate(it.product(outputFields, outputStats)):
    outName = '{0:s}_{1:s}.h5'.format(combos[0],combos[1])
    
    store = pd.HDFStore(outName) 

    for j in range(numTempBins):
        label = tempLabels[j]
        
        df = pd.DataFrame(arrs[i][j],columns=denseLabels,index=expns)
        print(df)
        store[label] = df

    store.close()
    
#
#dfvrFrac = pd.DataFrame(vrFraction,columns=labels,index=expns)
#dfvrStd = pd.DataFrame(vrStd,columns=labels,index=expns)
#dfrFrac = pd.DataFrame(rMean,columns=labels,index=expns)
#dfrStd = pd.DataFrame(rStd,columns=labels,index=expns)
#
#dfvrFrac.to_hdf('vradFraction_mean.h5','data',mode='w')
#dfvrStd.to_hdf('vradFraction_std.h5','data',mode='w')
#dfrFrac.to_hdf('r_mean.h5','data',mode='w')
#dfrStd.to_hdf('r_std.h5','data',mode='w')
#


        

