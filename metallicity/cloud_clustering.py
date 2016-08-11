

from __future__ import print_function
import numpy as np
import scipy.linalg as sl
import pandas as pd
import sklearn.cluster as skc
import sklearn.metrics as sm
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler


galNum = 25
a = '0.490'
numClouds = 3
loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5
baseLoc = '/home/jacob/research/velas/vela2b/vela{0:d}/a{1:s}/'
baseLoc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:s}/'

galNums = range(21,30)
expns = ['0.490']*len(galNums)
expns[galNums.index(24)] = '0.450'

for galNum, a in zip(galNums, expns):

    print('\nGalaxy = {0:d}:'.format(galNum))
    dataloc = baseLoc.format(galNum,a)
    boxfile = '{0:s}vela2b-{1:d}_GZa{2:s}.h5'.format(dataloc,galNum,a)
    # Read in the data and select out the cloud
    d = pd.read_hdf(boxfile, 'data')
    cloudInds = ( (d['temperature']<hiT) & (d['temperature']>loT) & 
                  (d['density']<hiN) & (d['density']>loN) )
    df = d[cloudInds]
    dataset = df[['x','y','z']]
    dloc = dataset.as_matrix()

    # Peform the kmean clustering
    #km = skc.KMeans(n_clusters=numClouds)
    #km.fit(dloc)
    #labels = km.labels_
    #results = pd.DataFrame([dataset.index,labels]).T

    # Check the fit with the Silhouette score
    #sm.silhouette_score(dloc, labels, metric='euclidean')

    # Compute the clustering using DBSCAN
    db = skc.DBSCAN(eps=0.3,min_samples=10).fit(dloc)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    print('Estimated number of clusters: %d' % n_clusters_)
    #print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(dloc, labels))



