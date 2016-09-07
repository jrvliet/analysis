
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

    binrange = [[0,200],[-6,-1]]

    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
    h, xedges, yedges, binnumber = st.binned_statistic_2d(locM['dist'],locM['logmetal'],
                                cloud['density'], statistic='count', bins=50,
                                range=binrange)

    print(np.sum(h))
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    mesh = ax1.pcolormesh(xedges, yedges, h)
    ax1.set_xlim(binrange[0])
    ax1.set_ylim(binrange[1])
    ax1.set_xlabel('Transverse Distance [kpc]')
    ax1.set_ylabel('Metallicity [MF]')
    cbar = plt.colorbar(mesh,ax=ax1,use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Counts',rotation=270,fontsize=12)
    
    h, xedges, yedges, binnumber = st.binned_statistic_2d(locM['dist'],locM['logmetal'],
                                cloud['density'], statistic='mean', bins=50, 
                                range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    mesh = ax2.pcolormesh(xedges, yedges, h)
    ax2.set_xlim(binrange[0])
    ax2.set_ylim(binrange[1])
    ax2.set_xlabel('Transverse Distance [kpc]')
    ax2.set_ylabel('Metallicity [MF]')
    cbar = plt.colorbar(mesh,ax=ax2,use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Mean Density',rotation=270,fontsize=12)
    

    sc = ax3.scatter(locM['dist'],locM['logmetal'], c=np.log10(cloud['density'].values), 
                    s=10, marker='o', alpha=0.01)
    binMeans, binEdges, binNumber = st.binned_statistic(locM['dist'], locM['logmetal'], 
                            statistic='mean', bins=20)
    binWidth = (binEdges[1]-binEdges[0])
    binCenters = binEdges[1:] - binWidth/2.
    ax3.hlines(binMeans, binEdges[:-1], binEdges[1:], colors='r', lw=2, label='Binned')
    cbar = plt.colorbar(sc, ax=ax3,use_gridspec=True)
    cbar.ax.set_ylabel('Density')
    ax3.set_xlim(binrange[0])
    ax3.set_ylim(binrange[1])
    ax3.set_xlabel('Distance from Line [kpc]')
    ax3.set_ylabel('Metallicity [MF]')
    


    fig.tight_layout()
    figName = 'vela2b-27_a{0:.3f}_Zgrandiant.png'.format(a)
    fig.savefig(figName, bbox_inches='tight', dpi=300)




















