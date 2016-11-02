
'''
Determines the net rotation of the cgm.
CGM is defined by all gas within 1.0 Rvir.
'''


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import sys

pd.options.mode.chained_assignment = None

# Temperature selections
tempEdges = np.linspace(3.5,6.5,4)
tempLabels = 'cool warm hot'.split()

# Density selections
denseEdges = np.linspace(-3.5,-6.5,7)
loNs = denseEdges[:-1]
hiNs = denseEdges[1:]
loNs = denseEdges[1:]
hiNs = denseEdges[:-1:]


dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.fullrot.h5'
rotfile = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

a = 0.490

rname = dataloc+rotfile.format(a)
with open(rname) as f:
    f.readline()
    rvir = float(f.readline().split()[3])

fname = dataloc+filename.format(a)
cgm = pd.read_hdf(fname,'data')

# Select out cgm
locIndex = (cgm['r']<1.0*rvir)
tempIndex = ((cgm['temperature']>10**tempEdges.min()) &
            (cgm['temperature']<10**tempEdges.max()))
densIndex = ((cgm['density']>10**denseEdges.min()) & 
            (cgm['density']<10**denseEdges.max()))
cgm = cgm[locIndex & tempIndex & densIndex]
print(cgm.shape)

cgm['speed'] = np.sqrt(cgm['x']**2 + cgm['y']**2 + cgm['z']**2)

basefields = ['x','y','z','vx','vy','vz','r','theta','phi','vr','vtheta','vphi']
coords = ['','Gal','Rot']

header = 'loN hiN'.split()
for coord in coords:
    for base in basefields:
        header.append('{0:s}{1:s}'.format(base,coord))


for i,temp in enumerate(tempLabels):
    print(temp)
    
    hdfFile = 'cgmRotation_{0:s}.h5'.format(temp)
    store = pd.HDFStore(hdfFile)

    loT = tempEdges[i]
    hiT = tempEdges[i+1]
    print(loT,hiT)

    tempIndex = (cgm['temperature']>10**loT) & (cgm['temperature']<10**hiT)

    mean = pd.DataFrame(index=range(len(loNs)),columns=header)
    stdev = pd.DataFrame(index=range(len(loNs)),columns=header)
    median = pd.DataFrame(index=range(len(loNs)),columns=header)
    minFrac = pd.DataFrame(index=range(len(loNs)),columns=header)
    maxFrac = pd.DataFrame(index=range(len(loNs)),columns=header)
    corr = pd.DataFrame(index=range(len(loNs)),columns=header)
    frames = [mean,stdev,median,minFrac,maxFrac,corr]
    
    for j,(loN,hiN) in enumerate(zip(loNs,hiNs)):
        print('\t{0:.1f} - {1:.1f}'.format(loN,hiN))

        for frame in frames:
            frame['loN'][j] = loN
            frame['hiN'][j] = hiN
            
        densIndex = (cgm['density']>10**loN) & (cgm['density']<10**hiN)
    
        gas = cgm[tempIndex & densIndex]
        print(gas.shape)

        for base in basefields:
            for coord in coords:
                field = '{0:s}{1:s}'.format(base,coord)
#                print('\t\t{0:s}'.format(field))
                
                mean[field][j] = gas[field].mean()
                stdev[field][j] = gas[field].std()
                median[field][j] = gas[field].median()
                minFrac[field][j] = (gas[field]/gas['speed']).min()
                maxFrac[field][j] = (gas[field]/gas['speed']).max()
                corr[field][j] = spearmanr(gas[field],gas['r'])[0]
        

    
    store['mean'] = mean
    store['stdev'] = stdev
    store['median'] = median
    store['minFrac'] = minFrac
    store['maxFrac'] = maxFrac
    store['spearman'] = corr
    
    store.close()
        




























