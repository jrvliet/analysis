
'''
This code is to investigate what happens to the warm gas in the filament during
the starburst (a=0.32-0.35). The mean metallicty of this gas somehow drops during
the starburst, which doesn't make any sense. This code will plot out the spatial
location of this gas with the hope that it will help 
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as st
import matplotlib.ticker as ticker

def fmt(x,pos):
    return '{0:.2f}'.format(x)

def mkHisto(ax,x,y,z,xlabel,ylabel,clabel='SNII'):

    numbins = 50
    stat = 'median'

    if 'z' in ylabel:
        binrange = [[-3,3],[0,6]]
    else:
        binrange = [[-3,3],[-3,3]]

    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,y,
                                statistic=stat,bins=numbins,range=binrange)
        
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)

    mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True,
                        format=ticker.FuncFormatter(fmt)) 
    



dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expns = np.arange(0.200,0.500,0.01)

loT,hiT = 4.5,5.5
lowestN,highestN = -5.5,-2.5
numDenseBins = 6
denseEdges = np.linspace(lowestN,highestN,numDenseBins+1)

for a in expns:
    print(a)

    fname = dataloc+filename.format(a)
    rname = dataloc+rotname.format(a)

    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])
    
    df = pd.read_hdf(fname,'data')

    spaceInd = (df['theta']<80)
    tempInd = (df['temperature']<10**hiT) & (df['temperature']>10**loT)

    fig,axes = plt.subplots(numDenseBins,3,figsize=(15,numDenseBins*5))
    
    for i in range(numDenseBins):

        loN = denseEdges[i]
        hiN = denseEdges[i+1]

        denseInd = (df['density']>10**loN) & (df['density']<10**hiN)

        cloud = df[spaceInd & tempInd & denseInd]
        
        # Plot out the three rotated coordinates plots
        ax = axes[i]
        
        mkHisto(ax[0],cloud['xRot']/rvir,cloud['yRot']/rvir,cloud['SNII'],'x','y')
        mkHisto(ax[1],cloud['xRot']/rvir,cloud['zRot']/rvir,cloud['SNII'],'x','z')
        mkHisto(ax[2],cloud['yRot']/rvir,cloud['zRot']/rvir,cloud['SNII'],'y','z')
        
        denseLabel = '{0:.1f}$<$log(nH)$<${1:.1f}\ny'.format(loN,hiN)
        ax[0].set_ylabel(denseLabel)

    fig.tight_layout()
    s = 'metallicity_drop_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)

         





