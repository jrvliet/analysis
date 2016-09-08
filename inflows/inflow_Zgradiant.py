
'''
Determines the profile of the metallicity gradiance 
of the ideal inflow in vela2b-27
'''

from __future__ import print_function
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.linalg as sl
import scipy.stats as st


def mkHist(ax,x,y,z,stat,binrange,xlab,ylab,zlab):
    
    h, xedges, yedges, binnumber = st.binned_statistic_2d(x,y,z, 
                                statistic=stat, bins=50, range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    mesh = ax.pcolormesh(xedges, yedges, h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(zlab,rotation=270,fontsize=12)
    
pd.options.mode.chained_assignment = None

galNum = 27
expns = [0.490]

# Define filenames
baseloc = '/home/jacob/research/velas/vela2b/vela{0:d}/a{1:.3f}/'
basename = 'vela2b-{0:d}_GZa{1:.3f}.h5'

# Define cloud parameters
loT, hiT = 10**3.25, 10**4.5
loN, hiN = 10**-6.25, 10**-2.25


for a in expns:
    
    loc = baseloc.format(galNum,a)
    fname = basename.format(galNum,a)

    df = pd.read_hdf(loc+fname, 'data')

    index = ( (df['temperature']<hiT) & (df['temperature']>loT) &
            (df['density']<hiN) & (df['density']>loN)  &
            (df['x']<0) & (df['z']>0) & (np.abs(df['y'])<300) )

    cloud = df[index]

    cloudLoc = cloud[['x', 'y', 'z']]
    locM = cloudLoc - cloudLoc.mean()
    covar = locM.cov()
    (u,s,v) = sl.svd(covar)
    
    u0 = u[0,0]
    u1 = u[1,0]
    u2 = u[2,0]
    
    # Get the distance from the line to each cell
    locM['a0'] = u2*locM['y'] - u1*locM['z']
    locM['a1'] = u0*locM['z'] - u2*locM['x']
    locM['a2'] = u1*locM['x'] - u0*locM['y']
    locM['dist'] = np.sqrt(locM['a0']**2 + locM['a1']**2 + locM['a2']**2)
    
    locM['metal'] = cloud['SNII'] + cloud['SNIa']
    locM['logmetal'] = np.log10(locM['metal'])
    cloud['r'] = np.sqrt( cloud['x']**2 + cloud['y']**2 + cloud['z']**2 )

    binrange = [[0,200],[-6,-1]]

    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))

    xlab = 'Transverse Distance [kpc]'
    ylab = 'Metallicity [MF]'
    mkHist(ax1,locM['dist'],locM['logmetal'],cloud['density'],'count',binrange,
            xlab,ylab,'Counts')
    mkHist(ax2,locM['dist'],locM['logmetal'],cloud['density'],'mean',binrange,
            xlab,ylab,'Mean Density')
    mkHist(ax3,locM['dist'],locM['logmetal'],cloud['r'],'mean',binrange,
            xlab,ylab,'Log (Mean Distance [kpc]) ')

    fig.tight_layout()
    figName = 'vela2b-27_a{0:.3f}_Zgrandiant.png'.format(a)
    fig.savefig(figName, bbox_inches='tight', dpi=300)




















