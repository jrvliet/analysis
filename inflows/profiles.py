
'''
Plots various properties as a function of distance from rotated z-axis
'''

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

pd.options.mode.chained_assignment = None

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

fields = 'density temperature SNII speed'.split()
markers = 's ^ d h o x'.split()

expns = np.arange(0.200,0.500,0.01)

numTempBins = 3
lowestT,highestT = 3.5,6.5
tempEdges = np.linspace(lowestT,highestT,numTempBins+1)
tempLabels = 'cool warm hot'.split()

lowestN,highestN = -5.5,-2.5
numDenseBins = 6
denseEdges = np.linspace(lowestN,highestN,numDenseBins+1)

denseLabels = []
for i in range(numDenseBins):
    denseLabels.append('{0:.1f}$<$n$<${1:.1f}'.format(denseEdges[i],denseEdges[i+1]))

for a in expns:

    print(a)

    fname = dataloc + filename.format(a)
    rname = dataloc + rotname.format(a)

    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    df = pd.read_hdf(fname,'data')

    df['rho'] = np.sqrt(df['xRot']**2 + df['yRot']**2)/rvir
    df['speed'] = np.sqrt(df['vxRot']**2 + df['vyRot']**2 + df['vzRot']**2)
    
    loT,hiT = tempEdges[0],tempEdges[1]
    tempLabel = tempLabels[0]
    tempInd = (df['temperature']<10**hiT) & (df['temperature']>10**loT)

    spaceInd = df['theta']<80
    
    fig,axes = plt.subplots(2,2,figsize=(10,10))
    
    for i in range(numDenseBins):

        loN = denseEdges[i]
        hiN = denseEdges[i+1]

        denseInd = (df['density']>10**loN) & (df['density']<10**hiN)

        cloud = df[denseInd & spaceInd & tempInd]
        
        for ax,field in zip(axes.flatten(),fields):

            h,binEdges,binnumber = st.binned_statistic(cloud['rho'],cloud[field],
                                    statistic='mean',bins=50,range=[0,4])
            binMid = []
            for j in range(len(h)):
                binMid.append((binEdges[j]+binEdges[j+1])/2.)
            ax.plot(binMid,h,marker=markers[i],label=denseLabels[i])
            ax.set_title(field)
    

    for ax in axes.flatten()[:-1]:
        ax.set_yscale('log')
    ax.legend(loc='best')
    fig.tight_layout()

    s = 'profile_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)
    













