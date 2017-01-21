
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


virialLoc = '/home/jacob/research/code/analysis/inflows/'
virialFile = 'vela2b-27_virialProperties.h5'
virial = pd.read_hdf(virialLoc+virialFile,'data')

loc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'

expns = np.arange(0.200,0.55,0.01)

lowestN,highestN = -5.5,-2.5
lowestT,highestT = 3.5,6.5
numNbins,numTbins = 3,3
nBins = np.linspace(lowestN,highestN,numNbins+1)
tBins = np.linspace(lowestT,highestT,numTbins+1)

tempLabels = 'cool warm hot'.split()
denseLabels = 'diffuse mid dense'.split()

for a in expns:

    print(a)

    rvir = virial['rvir'][np.isclose(virial['a'],a)].values
    print(rvir)
    fname = loc+filename.format(a)
    
    df = pd.read_hdf(fname,'data')

    spaceInds = df['theta']<80
    
    fig,axes = plt.subplots(3,3,figsize=(15,15))

    for i in range(numTbins):
        
        loT,hiT = tBins[i],tBins[i+1]
        tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
        
        for j in range(numNbins):

            loN,hiN = nBins[j],nBins[j+1]
            denseBins = (df['density']>10**loN) & (df['density']<10**hiN)

            fil = df[spaceInds & tempInds & denseBins]
            fil = fil.sample(frac=0.1)
    
            ax = axes[i,j]
            
            x = fil['yRot']/rvir
            y = fil['zRot']/rvir
            ax.scatter(x,y,marker='o',alpha=0.05)
            ax.set_xlim([-3,3])
            ax.set_ylim([0,4])
            
            tempLabel = tempLabels[i]
            denseLabel = denseLabels[j]
            ax.set_title('{0:s},{1:s}'.format(tempLabel,denseLabel))

    fig.tight_layout()
    s = 'filament_simple_a{0:d}.png'.format(int(a*1000))
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


