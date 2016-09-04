
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

pd.options.mode.chained_assignment = None

loT, hiT = 10**3.25, 10**4.5
loN, hiN = 10**-6.25, 10**-2.25
numbins = 20

dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/a0.490/'
dataloc = '/home/jacob/research/velas/vela2b/vela27/a0.490/'
fname = '{0:s}vela2b-27_GZa0.490.h5'.format(dataloc)

df = pd.read_hdf(fname, 'data')

index = ( (df['temperature']>loT) & (df['temperature']<hiT) & 
            (df['density']>loN) & (df['density']<hiN) &
            (df['x']>0) & (df['z']>0) & (np.abs(df['y'])<300) )

cloud = df[index]
cloudLoc = cloud[['x','y','z']]

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

# Determine the distance from each cell in the cloud
# to the this line
locM['a0'] = u2*locM['y'] - u1*locM['z']
locM['a1'] = u0*locM['z'] - u2*locM['x']
locM['a2'] = u1*locM['x'] - u0*locM['y']
locM['dist'] = np.sqrt(locM['a0']**2 + locM['a1']**2 + locM['a2'])

print(locM['dist'].min())
print(locM['dist'].max())
print(locM['dist'].mean())

out = pd.cut(locM['dist'], bins=numbins)
counts = pd.value_counts(out, sort=False)
print(counts)
print(type(counts))
counts.to_csv('counts.out',sep=' ',index=True,header=False)

(mu, sigma) = st.norm.fit(locM['dist'])
print('mu = {0:.3f}\nsigma = {1:.3f}'.format(mu,sigma))

# Plot historgram of dist
fig, ax = plt.subplots(figsize=(5,5))
n, bins, patches = ax.hist(locM['dist'], bins=numbins, range=[locM['dist'].min(),locM['dist'].max()],histtype='step',normed=True)
y = mlab.normpdf(bins, mu, sigma)
l = ax.plot(bins, y, 'r--', linewidth=2)

#pd.DataFrame.hist(locM,column='dist',ax=ax,bins=numbins,color='r')
ax.set_xlabel('Distance from line')
ax.set_title(r'$\mu = {0:.3f},\ \sigma = {1:.3f}$'.format(mu,sigma))

s = 'vela2b-27_inflow_coherence.png'
fig.savefig(s, bbox_inches='tight', dpi=300)




















