
# coding: utf-8

from __future__ import print_function
import numpy as np
import scipy.linalg as sl
import pandas as pd
import sklearn.cluster as skc
import sklearn.metrics as sm
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from collections import Counter
import itertools


galNum = 27
a = '0.490'
loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5
baseLoc = '/home/jacob/research/velas/vela2b/vela{0:d}/a{1:s}/'
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



# Compute the clustering using DBSCAN
eps = 0.5
min_samples=10
db = skc.DBSCAN(eps=eps,min_samples=min_samples).fit(dloc)
predict = skc.DBSCAN(eps=0.3,min_samples=10).fit_predict(dloc)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
n_clusters = len(set(labels)) - (1 if -1 in labels else 0)


# In[7]:

print('Number of cells in dloc = {0:d}'.format(len(dloc)))
print(' = {0:d}'.format(len(dloc)-sum(core_samples_mask)))
print('Number of clusters = {0:d}'.format(n_clusters))
Counter(labels)
u,c = np.unique(labels, return_counts=True)
for un,co in zip(u,c):
    print( un, co)


# In[ ]:

epsList = np.linspace(0.1,100.0,50)
msList = np.linspace(100,1000,10)
msList[0] = 1
with open('dbscan_params_test.out','w') as f:
    f.write('EPS\tMinSamples\tN_clusters\tPercent_Contatined\tN_Outliers\n')
    s = '{0:.2f}\t{1:.0f}\t\t{2:d}\t\t{3:.2%}\t\t\t{4:d}\n'
    for x in itertools.product(epsList,msList):
        db = skc.DBSCAN(eps=x[0],min_samples=x[1]).fit(dloc)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        frac = sum(core_samples_mask)/float(len(dloc))
        outliers = -1*sum(labels[labels==-1])
        f.write(s.format(x[0],x[1],n_clusters,frac,outliers))


# In[1]:

eps, minSamples, nClusters, nOutliers, percent = [],[],[],[],[]
with open('dbscan_params_test.out','r') as f:
    f.readline()
    for line in f:
        l = line.split()
        eps.append(float(l[0]))
        minSamples.append(float(l[1]))
        nClusters.append(float(l[2]))
        percent.append(float(l[3].strip('%')))
        nOutliers.append(float(l[4]))


# In[4]:

fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
axes = (ax1,ax2,ax3)
nClust = np.log10(nClusters)
nOut = np.log10(nOutliers)
colors = [nClust, nOut, percent]
titles = ['Number of Clusters','Number of Outliers','Percent Contained']

for i in range(len(axes)):
    ax = axes[i]
    a = ax.scatter(eps,minSamples,c=colors[i],marker='o',s=20,edgecolor=None)
    cbar = plt.colorbar(a,ax=ax, use_gridspec=True)
    ax.set_xlabel('EPS')
    ax.set_ylabel('Min Samples')
    ax.set_title(titles[i])
fig.tight_layout()    
plt.savefig('dbscan_parameters_results.png',  bbox_inches='tight', dpi=300)


# In[ ]:



