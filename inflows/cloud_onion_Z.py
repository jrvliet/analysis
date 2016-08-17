
'''
Plots the metallicity of gas as a funciton of distance
from the line of best fit along the the "cloud"
'''

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import scipy.linalg as sl
import pandas as pd
import sklearn.cluster as skc
import sklearn.metrics as sm
import sys

def plotting3d(cloud, ax3da, ax3db, ax3dc):

    elev, azim = 90, 0
    colors = ['red','blue','green']
    ax3da.scatter(cloud['x'],cloud['y'],cloud['z'], marker='o',c=colors[j], alpha=0.01)
    ax3da.scatter(x,y,z,marker='x',c='black',edgecolor='white')
    ax3da.view_init(elev=elev, azim=azim)
    ax3da.set_xlabel('x')
    ax3da.set_ylabel('y')
    ax3da.set_zlabel('z')
    ax3da.set_title('Elev={0:d}, Azim={1:d}'.format(elev,azim))
    
    elev, azim = 0, 0
    ax3db.scatter(cloud['x'],cloud['y'],cloud['z'], marker='o',c=colors[j], alpha=0.01)
    ax3db.scatter(x,y,z,marker='x',c='black',edgecolor='white')
    ax3db.view_init(elev=elev, azim=azim)
    ax3db.set_xlabel('x')
    ax3db.set_ylabel('y')
    ax3db.set_zlabel('z')
    ax3db.set_title('Elev={0:d}, Azim={1:d}'.format(elev,azim))
    
    elev, azim = 90, 90
    ax3dc.scatter(cloud['x'],cloud['y'],cloud['z'], marker='o',c=colors[j], alpha=0.01)
    ax3dc.scatter(x,y,z,marker='x',c='black',edgecolor='white')
    ax3dc.view_init(elev=elev, azim=azim)
    ax3dc.set_xlabel('x')
    ax3dc.set_ylabel('y')
    ax3dc.set_zlabel('z')
    ax3dc.set_title('Elev={0:d}, Azim={1:d}'.format(elev,azim))


def plotting2d(ax, axb, locM, j):
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


plotting = 0

# Define cloud parameters
numClouds = 3
loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5

galNums = range(21,30)
expns = ['0.490']*len(galNums)
expns[galNums.index(24)] = '0.450'

baseLoc = '/home/jacob/research/velas/vela2b/vela{0:d}/a{1:s}/'
baseLoc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:s}/'

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
    silhouette = sm.silhouette_score(dloc, labels, metric='euclidean')

    print('\tClustering done')

    fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6))) = plt.subplots(3,2,figsize=(10,15))
    axes = [ax1,ax2,ax3,ax4,ax5,ax6]
    figb, axb = plt.subplots(1,1,figsize=(10,10))
    fig3d = plt.figure(figsize=(15,10))
    ax3da = fig3d.add_subplot(131,projection='3d')
    ax3db = fig3d.add_subplot(132,projection='3d')
    ax3dc = fig3d.add_subplot(133,projection='3d')

    # Loop through the seperate 'clouds'
    for j in range(numClouds):
        print('\t\tCloud {0:d}'.format(j+1))

        cloud = df[(results[1]==j).values]
        cLoc = cloud[['x','y','z']]
        
        # Determine the mean location of this cloud
        x = cLoc['x'].mean()
        y = cLoc['y'].mean()
        z = cLoc['z'].mean()

        # Plot the 3D 
        if plotting==1:
            plotting3d(cloud, ax3da, ax3db, ax3dc)
        
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

        if plotting==1:
            ax = axes[j]
            plotting2d(ax, axb, locM, j)
        
    if plotting==1:
        axb.legend()
        fig.tight_layout()
        fig3d.tight_layout()
        s1 = 'vela2b-{0:d}_cloud_onion.png'.format(galNum)
        s2 = s1.replace('onion','onion_all')
        s3d = s1.replace('onion','clustering')
        fig.savefig(s1, bbox_inches='tight',dpi=300)
        figb.savefig(s2, bbox_inches='tight',dpi=300)
        fig3d.savefig(s3d, bbox_inches='tight',dpi=300)

        plt.cla()
        plt.clf()

hdf5file = 'vela2b_cloud_loc_stats.h5'
df = pd.DataFrame(data,columns=header)
df.to_hdf(hdf5file,'data',mode='w')





