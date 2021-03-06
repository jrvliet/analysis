
# coding: utf-8

# In[10]:

from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
pd.options.mode.chained_assignment = None
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import sys


def fmt(x,pos):
    return '{0:.1f}'.format(x)

# In[16]:

def mkHist(ax,x,y,z,stat,xlabel,ylabel,lox,hix,loy,hiy):
    numbins = 50
    #stat = 'count'
    binrange = [[lox,hix],[loy,hiy]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,z,
                              statistic=stat, bins=numbins,
                                 range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0

    # Convert the cells from total mass to fraction of total masss of all cells
    totalMass = float(h.sum())
    h = h/totalMass

    vmax = -1
    vmin = -6
    #if vmax<h.max():
    #    print('Exceeds\t',h.max())
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    print(h.max(),np.percentile(h,10))

    #mesh = ax.pcolormesh(xedges,yedges,h,vmin=0,vmax=vmax)
    mesh = ax.pcolormesh(xedges,yedges,h,vmin=vmin,vmax=vmax,cmap='jet')
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

loc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
loc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'
expns = np.arange(0.210,0.540,0.01)


sfrLoc = '/home/jacob/research/code/analysis/sfr/'
sfrFile = 'vela2b-27_sfr.csv'
sfr = pd.read_csv(sfrLoc+sfrFile)

# In[4]:

loT,hiT = 3.5,6.5
loN,hiN = -5.5,-2.5

header = 'a nMax tMax'.split()
maxMetal = np.zeros((len(expns),len(header)))
maxMetal = pd.DataFrame(maxMetal,columns=header)

print(sfr)
# In[18]:

for i,a in enumerate(expns):
    print(a)

    rname = loc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = loc+filename.format(a)
    df = pd.read_hdf(fname,'data')
    df['r'] = np.sqrt(df['xRot']**2 + df['yRot']**2 + df['zRot']**2)
    tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    densInds = (df['density']>10**loN) & (df['density']<10**hiN)
    spacInds = (df['theta']<80)
    fil = df[tempInds & densInds & spacInds]
    
    fil['mass'] = cellMass(fil['density'],fil['cell_size'])
    fil['metalMass'] = fil['mass']*(fil['SNII']+fil['SNIa'])
    





    #ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
    #ax2 = plt.subplot2grid((3,1), (2,0))

    fig = plt.figure(figsize=(5,7))
    gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
    gs.update(hspace=0.25)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    h,xedges,yedges,binnumber = mkHist(ax1,np.log10(fil['density']),
                                np.log10(fil['temperature']),fil['metalMass'],
                                'sum','Density','Temperature',loN,hiN,loT,hiT)
    ax1.set_title('{0:.3f}'.format(a))

    ax2.plot(sfr['a'],sfr['sfr'],color='black')
    ax2.set_xlabel('Expansion Parameter')
    ax2.set_ylabel(r'SFR [M$_{\odor}$ yr$^{-1}$]')
    ymin,ymax = ax2.get_ylim()
    ax2.scatter(sfr['a'].ix[i],sfr['sfr'].ix[i],color='red',marker='s')
    
    s = 'metalMassHist_a{0:.3f}_noCGM_fraction.png'.format(a)
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

    
hdfFile = 'maxMetalPhase_noCGM_fraction.h5'
maxMetal.to_hdf(hdfFile,'data',mode='w')



