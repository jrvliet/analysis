
'''
Plots the metallicity of gas as a funciton of distance
from the line of best fit along the the "cloud"
'''

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sl
import pandas as pd
import sklearn.cluster as skc
import sys

# Define cloud parameters
numClouds = 3
loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5

galNums = range(21,30)
expns = ['0.490']*len(galNums)
expns[galNums.index(24)] = '0.450'

baseLoc = '/home/jacob/research/velas/vela2b/vela{0:d}/a{1:s}/'

header = ['galNum','x1','y1','z1','r1','theta1','phi1','sn1','vrad1',
        'x2','y2','z2','r2','theta2','phi2','sn2','vrad2'
        'x3','y3','z3','r3','theta3','phi3','sn3','vrad3']

data = np.zeros((len(galNums),len(header)))

for i,(galNum, a) in enumerate(zip(galNums, expns)):

    print('vela2b-{0:d}'.format(galNum)) 
    data[i,0] = galNum

    # Read in data
    dataloc = baseLoc.format(galNum,a)
    boxfile = '{0:s}vela2b-{1:d}_GZa{2:s}.h5'.format(dataloc,galNum,a)

    d = pd.read_hdf(boxfile, 'data')
    print('\tData read in')

    # Select the cells in the cloud
    cloudInds = ( (d['temperature']<hiT) & (d['temperature']>loT) & 
                  (d['density']<hiN) & (d['density']>loN) )
    df = d[cloudInds]
    dataset = df[['x','y','z']]
    dloc = dataset.as_matrix()
    print('\tCloud selected')

    # Sort the cloud cells into 3 groups using k-means clustering with n=3
    km = skc.KMeans(n_clusters=numClouds)
    km.fit(dloc)
    labels = km.labels_
    results = pd.DataFrame([dataset.index,labels]).T
    print('\tClustering done')

    fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(5,15))
    axes = [ax1,ax2,ax3]
    figb, axb = plt.subplots(1,1,figsize=(10,10))
    # Loop through the seperate 'clouds'
    for j in range(numClouds):
        print('\t\tCloud {0:d}'.format(j+1))

        cloud = df[(results[1]==j).values]
        cLoc = cloud[['x','y','z']]
        
        # Determine the mean location of this cloud
        x = cLoc['x'].mean()
        y = cLoc['y'].mean()
        z = cLoc['z'].mean()

        # Convert to spherical coordinates
        cloud.loc[:,'dist'] = np.sqrt( cloud['x']**2 + cloud['y']**2 + cloud['z']**2)
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arctan2(y,x)*(180.0/np.pi)
        phi = np.arccos(z/r)*(180.0/np.pi)

        # Calculate the mean radial velocity
        cloud.loc[:,'radvel'] = (cloud['x']*cloud['vx'] + cloud['y']*cloud['vy'] + 
                            cloud['z']*cloud['vz'] ) / cloud['dist']
        
        row = j*6
        data[i,row+1] = x
        data[i,row+2] = y
        data[i,row+3] = z
        data[i,row+4] = r
        data[i,row+5] = theta
        data[i,row+6] = phi
        data[i,row+7] = cloud['SNII'].mean()
        data[i,row+8] = cloud['radvel'].mean()
        
        # Normalize the cells by subtracting out their mean
        locM = cLoc - cLoc.mean()

        # Get the covariance matrix
        covar = locM.cov()

        # Determine the single value decomposition of the covariance matrix
        (u,s,v) = sl.svd(covar)

        # The first column in u is the directional vector of the line
        # of best fit for the cloud
        u0 = u[0,0]
        u1 = u[1,0]
        u2 = u[2,0]

        # To find the distance from each cell to this line, use:
        # dist = | (v2-v1) x (v1-v0) | 
        #        ---------------------
        #             | v2 - v1 |
        # where:
        # v0 = vector to the cell
        # v1 = vector to a point at the start of the line
        # v2 = vector to a point at the end of the line
        #
        # Here, the coordinates are the centered coordinates, so
        # v1 = 0 and v2 is given by <u0,u1,u2> and has length 1
        # The coordinates of the cells need to be pulled from locM, not cloud
        # With these considerations, the formula simplifies to:
        # dist = | v2 x (-v0) |
        # Define:
        # a = v2 x (-v0)
        # where:
        # a0 = u2*y - u1*z
        # a1 = u0*z - u2*x
        # a2 = u1*x - u0*y
        # dist = sqrt( a0**2 + a1**2 + a2**2 )
        # Add columns to the locM dataframe with these properties
        locM['a0'] = u2*locM['y'] - u1*locM['z']
        locM['a1'] = u0*locM['z'] - u2*locM['x']
        locM['a2'] = u1*locM['x'] - u0*locM['y']
        locM['dist'] = np.sqrt(locM['a0']**2 + locM['a1']**2 + locM['a2']**2)

        # Get the metallicity of each of these cells
        locM['metal'] = cloud['SNII'] + cloud['SNIa']

        ax = axes[j]
        ax.plot(locM['dist'], locM['metal'], '.', alpha=0.01)
        ax.set_xlabel('Distance [kpc]')
        ax.set_ylabel('Metals [mass fraction]')
        ax.set_yscale('log')
        ax.set_title('Cloud {0:d}'.format(j+1))    

        axb.plot(locM['dist'], locM['metal'], '.', alpha=0.01, label='Cloud {0:d}'.format(j+1))
        axb.set_xlabel('Distance [kpc]')
        axb.set_ylabel('Metals [mass fraction]')
        axb.set_yscale('log')
        axb.set_title('Cloud {0:d}'.format(j+1))    

    axb.legend()
    fig.tight_layout()
    s1 = 'vela2b-{0:d}_cloud_onion.png'.format(galNum)
    s2 = s1.replace('onion','onion_all')
    fig.savefig(s1, bbox_inches='tight',dpi=300)
    figb.savefig(s2, bbox_inches='tight',dpi=300)

hdf5file = 'vela2b_cloud_loc_stats.h5'
df = pd.DataFrame(data,columns=header)
df.to_hdf(hdf5file,'data',mode='w')





