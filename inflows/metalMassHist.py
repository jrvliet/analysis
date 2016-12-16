
# coding: utf-8

# In[10]:

from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
pd.options.mode.chained_assignment = None
import matplotlib.ticker as ticker


def fmt(x,pos):
    return '{0:.1f}'.format(x)

# In[16]:

def mkHist(ax,x,y,z,stat,xlabel,ylabel,lox,hix,loy,hiy):
    numbins = 500
    #stat = 'count'
    binrange = [[lox,hix],[loy,hiy]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,z,
                              statistic=stat, bins=numbins,
                                 range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)

    mesh = ax.pcolormesh(xedges,yedges,h,vmin=0)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True,
                       format=ticker.FuncFormatter(fmt) )
    return h,xedges,yedges,binnumber


# In[5]:

def cellMass(n,l):
    kpc2km = 3.086e16 
    pc2cm = 3.086e18 
    mH = 1.6737e-24 
    mSun = 1.989e33 
    s2yr = 3.154e7 
    return n*mH*(l*pc2cm)**3 / mSun


# In[19]:

loc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
expns = np.arange(0.200,0.550,0.01)


# In[4]:

loT,hiT = 3.5,6.5
loN,hiN = -5.5,-2.5

header = 'a nMax tMax'.split()
maxMetal = np.zeros((len(expns),len(header)))
maxMetal = pd.DataFrame(maxMetal,columns=header)

# In[18]:

for i,a in enumerate(expns):
    print(a)
    fname = loc+filename.format(a)
    df = pd.read_hdf(fname,'data')
    tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    densInds = (df['density']>10**loN) & (df['density']<10**hiN)
    spacInds = (df['theta']<80)
    fil = df[tempInds & densInds & spacInds]
    
    fil['mass'] = cellMass(fil['density'],fil['cell_size'])
    fil['metalMass'] = fil['mass']*(fil['SNII']+fil['SNIa'])
    
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    h,xedges,yedges,binnumber = mkHist(ax,np.log10(fil['density']),
                                np.log10(fil['temperature']),fil['metalMass'],
                                'sum','Density','Temperature',loN,hiN,loT,hiT)
    ax.set_title('{0:.3f}'.format(a))
    s = 'metalMassHist_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)

    # Get the location of the maximum metal mass
    tMaxInd,nMaxInd = np.unravel_index(h.argmax(), h.shape)
    tMax = yedges[tMaxInd]
    nMax = xedges[nMaxInd]

    # Record these values
    maxMetal['a'].ix[i] = a
    maxMetal['nMax'].ix[i] = nMax
    maxMetal['tMax'].ix[i] = tMax

    
hdfFile = 'maxMetalPhase.h5'
maxMetal.to_hdf(hdfFile,'data',mode='w')



