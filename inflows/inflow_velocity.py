

'''
Determines the coherence of the inflows
Focus on the poster child inflow in vela2b-27
at a=0.490
'''


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.linalg as sl
import scipy.stats as st
import sys

def mkHist(ax,x,lims,numbins,label,xlab, log=False):

    ax.hist(x,bins=numbins,range=lims,histtype='step',log=log,label=label)
    ax.set_xlabel(xlab)
    ax.set_ylabel('Counts')

pd.options.mode.chained_assignment = None

loT, hiT = 10**3.25, 10**4.5
loN, hiN = 10**-6.25, 10**-2.25
numbins = 200

expns = [0.490]

baseloc = '/home/jacob/research/velas/vela2b/vela27/a{0:.3f}/'
basename = '{0:s}vela2b-27_GZa{1:.3f}.h5'

for a in expns:

    print('a = {0:.3f}'.format(a))
    
    loc = baseloc.format(a)
    fname = basename.format(loc,a)
    
    df = pd.read_hdf(fname, 'data')

    index = ( (df['temperature']>loT) & (df['temperature']<hiT) & 
                (df['density']>loN) & (df['density']<hiN) &
                (df['x']>0) & (df['z']>0) & (np.abs(df['y'])<300) )


    # Start working with velocities
    df['speed'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2 )

    # Determine the speed of each cell
    loSpeed = df['speed'].min()
    hiSpeed = df['speed'].max()

    # Calculate the spherical coordinates
    df['r'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2 )
    df['theta'] = np.arctan2(df['y'],df['x'])
    df['phi'] = np.arccos(df['z']/df['r'])
    
    # Calculate the spherical velocities
    df['vr'] = (df['x']*df['vx'] + df['y']*df['vy'] + 
                    df['z']*df['vz']) / df['r']
    df['vtheta'] = ( (df['x']*df['vy'] - df['y']*df['vx']) /
                        (df['x']**2 + df['y']**2 ) )
    df['vphi'] = ( (-1 / np.sqrt( 1 - (df['z']/df['r'])**2)) * 
                        (df['vz']/df['r'] - df['z']*df['vr']/df['r']**2) )

    cloud = df[index]
    cloudLoc = cloud[['x','y','z']]
    print('Number of samples = {0:d}'.format(len(cloud)))

    # Perform a PCA to find the line of best fit
    # Start by normalizing the cells
    locM = cloudLoc - cloudLoc.mean()

    # Get the covariance matrix
    covar = locM.cov()

    # Determine the single value decomposition of the covariance matrix
    (u,s,v) = sl.svd(covar)

    # The first column of u is the directional vector of the line
    # of best fit for the cloud
    u0 = u[0,0]
    u1 = u[1,0]
    u2 = u[2,0]

    adotb = cloud['vx']*u0 + cloud['vy']*u1 + cloud['vz']*u2
    bdotb = u0**2 + u1**2 + u2**2
    factor = adotb/bdotb

    cloud['along0'] = factor*u0
    cloud['along1'] = factor*u1
    cloud['along2'] = factor*u2
    cloud['along'] = np.sqrt(cloud['along0']**2 + 
                        cloud['along1']**2 + cloud['along2']**2)

    cloud['perp0'] = cloud['vx'] - cloud['along0']
    cloud['perp1'] = cloud['vy'] - cloud['along1']
    cloud['perp2'] = cloud['vz'] - cloud['along2']
    cloud['perp'] = np.sqrt(cloud['perp0']**2 + 
                        cloud['perp1']**2 + cloud['perp2']**2)


    # Plot historgram of dist
    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(15,15))
    log = False

    # Plot cartesian velocities
    lims = [cloud['vx'].min(), cloud['vx'].max()]
    mkHist(ax1, cloud['vx'], lims, numbins, 'cloud', 'vx', log)

    lims = [cloud['vy'].min(), cloud['vy'].max()]
    mkHist(ax2, cloud['vy'], lims, numbins, 'cloud', 'vy', log)

    lims = [cloud['vz'].min(), cloud['vz'].max()]
    mkHist(ax3, cloud['vz'], lims, numbins, 'cloud', 'vz', log)

    # Plot spherical velocities
    log = True
    lims = [cloud['vr'].min(), cloud['vr'].max()]
    mkHist(ax4, cloud['vr'], lims, numbins, 'cloud', 'vr', log)
    mkHist(ax4, df['vr'], lims, numbins, 'all', 'vr', log)

    lims = [cloud['vphi'].min(), cloud['vphi'].max()]
    mkHist(ax5, cloud['vphi'], lims, numbins, 'cloud', 'vphi', log)
    mkHist(ax5, df['vphi'], lims, numbins, 'all', 'vphi', log)

    lims = [cloud['vtheta'].min(), cloud['vtheta'].max()]
    mkHist(ax6, cloud['vtheta'], lims, numbins, 'cloud', 'vtheta', log)
    mkHist(ax6, df['vtheta'], lims, numbins, 'all', 'vtheta', log)
    
    ax4.legend()
    ax5.legend()
    ax5.legend()
    
    # Plot the speeds
    log = False
    lims = [cloud['speed'].min(), cloud['speed'].max()]
    mkHist(ax7, cloud['speed'], lims, numbins, 'cloud', 'Speed')

    lims = [cloud['along'].min(), cloud['along'].max()]
    mkHist(ax8, cloud['along'], lims, numbins, 'cloud', 'Along')

    lims = [cloud['perp'].min(), cloud['perp'].max()]
    mkHist(ax9, cloud['perp'], lims, numbins, 'cloud', 'V Perp')
    
    fig.tight_layout()
    s = 'vela2b-27_a{0:.3f}_inflowVelocity.png'.format(a)
    fig.savefig(s, bbox_inches='tight', dpi=300)







