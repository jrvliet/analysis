
'''
Determines the profile of metallicity in possible 
inflows
'''


from __future__ import print_function
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import sklearn.cluster as skc
import scipy.linalg as sl
import scipy.stats as st

pd.options.mode.chained_assignment = None  # default='warn'
colors = ['blue','red','green','black','cyan','coral','purple','goldenrod','gray','orange']
markers = ['o','^','d','s','v']

# Define which galaxies to use
galNums = range(21,30)
expns = ['0.490']*len(galNums)
expns[galNums.index(24)] = '0.450'

# Define cloud parameters
loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5

# Define DBSCAN parameters
eps = 60        # Maximum distance between cells in a cluster
minPart = 500   # Minimum number of cells to be called a cluster

# Define filenames
baseloc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:s}/'
filename = 'vela2b-{0:d}_GZa{1:s}.h5'

outfile = 'vela2b_onion_DBSCAN_eps{0:d}_min{1:d}_results.dat'.format(eps,minPart)
f = open(outfile, 'w')
header = 'GalNum\tnCluster\tnOutliers\tnTotal\t%In\n'
f.write(header)
ss = '{0:d}\t{1:d}\t\t{2:d}\t\t{3:d}\t{4:.1%}\n'

fitfile = 'vela2b_onion_DBSCAN_eps{0:d}_min{1:d}_fit.dat'.format(eps,minPart)
fitf = open(fitfile, 'w')
header = 'GalNum\tCluster\tSlope\tInter\tR\tP\tStd_err\n'
fitf.write(header)
fits = '{0:d}\t{1:d}\t{2:.3e}\t{3:.3e}\t{4:.3f}\t{5:.3f}\t{6:.3f}\n'

velfile = 'vela2b_onion_DBSCAN_eps{0:d}_min{1:d}_velocity.h5'.format(eps,minPart)
velheader = ['GalNum','Cluster','vr_mean','vr_min','vr_max','vr_std','r_mean','r_min','r_max','r_std']
numcols = len(velheader)
veldf = np.zeros(numcols)


# Loop over galaxies
for galNum, a in zip(galNums, expns):

    dataloc = baseloc.format(galNum,a)
    fname = filename.format(galNum,a)

    # Read in the data
    df = pd.read_hdf(dataloc+fname, 'data')

    # Select out the cloud
    cloudInds = ( (df['temperature']<hiT) & (df['temperature']>loT) &
                  (df['density']<hiN) & (df['density']>loN) )
    cloud = df[cloudInds]

    # Select out only the coordinates for the fitting
    dloc = cloud[['x','y','z']].as_matrix()

    # Perform the DBSCAN
    db = skc.DBSCAN(eps=eps, min_samples=minPart).fit(dloc)

    # Pull out the results of the clustering analysis
    core_samples_mask = np.zeros_like(db.labels_,dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    nClusters = len(set(labels)) - (1 if -1 in labels else 0)

    # Write results to output
    nOutliers = list(labels).count(-1)
    nTotal = len(df)
    fraction = (nTotal-nOutliers)/float(nTotal)
    f.write(ss.format(galNum,nClusters,nOutliers,nTotal,fraction)) 

    # Print out to screen
    print('Galaxy Number = {0:d}\tNumber of Clusters = {1:d}'.format(galNum,nClusters))
    
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111, projection='3d')
    # Loop over each cluster
    for i in range(nClusters):
        
        ind = (labels==i)
        cluster = cloud[ind]
        clusterLoc = cluster[['x','y','z']]

        m = i%len(markers)
        c = i%len(colors)
        ax.scatter(cluster['x'],cluster['y'],cluster['z'],
                    marker=markers[m],color=colors[c],alpha=0.01)

        # Determine the radial properties of this cluster
        cluster['r'] = np.sqrt(cluster['x']**2 + cluster['y']**2 + cluster['z']**2)
        cluster['vr'] = (cluster['vx']*cluster['x'] + cluster['vy']*cluster['y'] + 
                            cluster['vz']*cluster['z']) / cluster['r']
        
        # Determine the mean parameters
        xmean = cluster['x'].mean()
        ymean = cluster['y'].mean()
        zmean = cluster['z'].mean()

        clust = np.zeros(numcols)
        clust[0] = galNum
        clust[1] = i
        clust[2] = cluster['vr'].mean()
        clust[3] = cluster['vr'].min()
        clust[4] = cluster['vr'].max()
        clust[5] = cluster['vr'].std()
        clust[6] = cluster['r'].mean()
        clust[7] = cluster['r'].min()
        clust[8] = cluster['r'].max()
        clust[9] = cluster['r'].std()

        veldf = np.vstack((veldf,clust))

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

        # Fit a line to the data
        slope, intercept, r_value, p_value, std_err = st.linregress(locM['dist'],locM['metal'])
        
        # Write the results to file
        fitf.write(fits.format(galNum,i,slope,intercept,r_value,p_value,std_err))

    # Save the plot to file
    ax.view_init(elev=0,azim=0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim([-300,300])
    ax.set_ylim([-300,300])
    ax.set_zlim([-300,300])
    fig.savefig('vela2b-{0:d}_eps{1:d}_min{2:d}_dbscanClusters.png'.format(galNum,eps,minPart),bbox_inches='tight',dpi=300)

f.close()
fitf.close()
df = pd.DataFrame(veldf,columns=velheader)
df.to_hdf(velfile,'data',mode='w')






