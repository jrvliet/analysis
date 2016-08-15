
'''
Determines the profile of metallicity in possible 
inflows
'''


from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.cluster as skc
import scipy.linalg as sl

# Define which galaxies to use
galNums = range(21,30)
expns = ['0.490']*len(galNums)
expns[galNums.index(24)] = '0.450'

# Define cloud parameters
loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5

# Define DBSCAN parameters
eps = 10        # Maximum distance between cells in a cluster
minPart = 100   # Minimum number of cells to be called a cluster

# Define filenames
baseloc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:s}/'
filename = 'vela2b-{0:d}_GZa{1:s}.h5'
outfile = 'vela2b_onion_DBSCAN_results.dat'
f = open(outfile, 'w')
header = 'GalNum\tnCluster\tnOutliers\nnTotal\t%In\n'
f.write(header)
s = '{0:d}\t{1:d}\t{2:d}\t{3:d}\t{4:.1%}\n'


# Loop over galaxies
for galNum, a in zip(galNums, expns):

    dataloc = baseloc.format(galNum,a)
    fname = filename.format(galNum,a)

    # Read in the data
    df = pd.read_hdf(dataloc+fname, 'data')

    # Select out the cloud
    cloudInds = ( (d['temperature']<hiT) & (d['temperature']>loT) &
                  (d['density']<hiN) & (d['density']>loN) )
    cloud = df[cloudInds]

    # Select out only the coordinates for the fitting
    dloc = cloud[['x','y','z']].as_matrix()

    # Perform the DBSCAN
    db = skc.DBSCAN(eps=eps, min_samples=minPart).fit(dloc)

    # Pull out the results of the clustering analysis
    core_samples_mask = np.zeros_like(db.labels_,dtype=bool)
    core_samples_mask[db.core_samples_indicies_] = True
    labels = db.labels_
    nClusters = len(set(labels)) - (1 if -1 in labels else 0)

    # Write results to output
    nOutliers = labels.count(-1)
    nTotal = len(df)
    fraction = (nTotal-nOutliers)/float(nTotal)
    f.write(s.format(galNun,nClusters,nOutliers,nTotal,fraction)) 
    
    # Loop over each cluster
    for i in range(nClusters):
        
        ind = (labels==i)
        cluster = cloud[ind]
        clusterLoc = cluster[['x','y','z']]

        # Determine the radial properties of this cluster
        cluster['r'] = np.sqrt(cluster['x']**2 + cluster['y']**2 + cluster['z']**2)
        cluster['vr'] = (cluster['vx']*cluster['x'] + cluster['vy']*cluster['y'] + 
                            cluster['vz']*cluster['z']) / cluster['r']
        
        # Determine the mean parameters
        xmean = cluster['x'].mean()
        ymean = cluster['y'].mean()
        zmean = cluster['z'].mean()
        rmean = cluster['r'].mean()
        vrmean = cluster['vr'].mean()

        # Perform PCA to find the line of best fit
        # Start by normalizing the cells
        locM = clusterLoc - clusterLoc.mean()

        # Get the covariance matrix
        covar = locM.cov()

        # Determine the single value decomposition of the covariance matrix
        (u,s,v) = sl.svd(covar)

        # The first column of u is the directional vector of the line
        # of best fit for the cloud
        u0 = u[0,0]
        u1 = u[1,0]
        u2 = u[2,0]
        
        # Determine the distance from each cell to this line
        # See comments in cloud_onion_Z.py for logic
        locM['a0'] = u2*locM['y'] - u1*locM['z']
        locM['a1'] = u0*locM['z'] - u2*locM['x']
        locM['a2'] = u1*locM['x'] - u0*locM['y']
        locM['dist'] = np.sqrt(locM['a0']**2 + locM['a1']**2 + locM['a2'])

        # Get the metallicity of each cell
        locM['metal'] = cluster['SNII'] + cluster['SNIa']


f.close()






