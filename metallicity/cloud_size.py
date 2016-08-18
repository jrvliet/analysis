
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import gc
import statsmodels as sm


dataloc = '/home/jacob/research/velas/vela2b/vela27/a0.490/'
#dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/a0.490/'
boxfile = '{0:s}vela2b-27_GZa0.490.h5'.format(dataloc)

d = pd.read_hdf(boxfile, 'data')

dense = np.arange(-5,-2,0.25)
loT, hiT = 10**4, 10**4.5
loN = 10**-5

baseInds = ( (d['temperature']<=hiT) & (d['temperature']>=loT) & 
              (d['x']<0) & (d['z']>0) & (np.abs(d['y'])<300))

fig = plt.figure(figsize=(20,15))
for i,n in enumerate(dense):
    print i
    if i>=1:
        loN = 10**(dense[i-1])
    hiN = 10**n
    cloudInds = (baseInds & (d['density']<=hiN) & (d['density']>=loN))
    print 'log(loN) = {0:.2f}\tlog(hiN) = {1:.2f}\tsum(cloudInds) = {2:d}'.format(np.log10(loN),np.log10(hiN),sum(cloudInds))
    cloud = d[cloudInds]
    ax = fig.add_subplot(3, 4, i+1, projection='3d')
    ax.scatter(cloud['x'], cloud['y'], cloud['z'], c='b', marker='o', alpha=0.01)
    ax.view_init(elev=90, azim=0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim([-350,0])
    ax.set_ylim([-150,150])
    ax.set_zlim([-50,350])
    ax.view_init(elev=30, azim=0)
    ax.set_title('{0:.2f} $<$ log(nH) $<$ {1:.2f}'.format(np.log10(loN),np.log10(hiN)))


fig.tight_layout()
plt.savefig('vela2b-27_a0.490_cloudSize.png', bbox_inches='tight', pdi=300)



