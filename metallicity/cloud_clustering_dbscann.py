

from __future__ import print_function
import numpy as np
import scipy.linalg as sl
import pandas as pd
import sklearn.cluster as skc
import sklearn.metrics as sm
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler


galNum = 27
a = '0.490'
loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5
baseLoc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:s}/'

dataloc = baseLoc.format(galNum,a)
boxfile = '{0:s}vela2b-{1:d}_GZa{2:s}.h5'.format(dataloc,galNum,a)

# Read in the data and select out the cloud
d = pd.read_hdf(boxfile, 'data')
cloudInds = ( (d['temperature']<hiT) & (d['temperature']>loT) & 
              (d['density']<hiN) & (d['density']>loN) )
df = d[cloudInds]
dataset = df[['x','y','z']]
dloc = dataset.as_matrix()

eps = 0.3
min_samples = range(2,11)

fig = plt.figure(figsize=(15,15))
for i in range(len(min_samples)):

    # Compute the clustering using DBSCAN
    db = skc.DBSCAN(eps=0.3,min_samples=10).fit(dloc)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    reults = 

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    ax = fig.add_subplot(3,3,i+1,projection='3d')





