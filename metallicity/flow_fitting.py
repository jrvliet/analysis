import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


dataloc = '/home/jacob/research/velas/vela2b/vela27/a0.490/'
dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/a0.490/'
boxfile = '{0:s}vela2b-27_GZa0.490.h5'.format(dataloc)


d = pd.read_hdf(boxfile, 'data')
loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5
print len(d)

cloudInds = ( (d['temperature']<hiT) & (d['temperature']>loT) & 
              (d['density']<hiN) & (d['density']>loN) &
              (d['x']<0) & (d['z']>0) & (np.abs(d['y'])<300))

cloud = d[cloudInds]

loc = cloud[['x','y','z']]
locMat = loc.as_matrix()
datamean = loc.mean(axis=0).as_matrix()

uu, dd, vv = np.linalg.svd(locMat, full_matrices=True)

print uu
print dd
print vv



