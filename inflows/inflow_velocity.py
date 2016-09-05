
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
import sys

pd.options.mode.chained_assignment = None

loT, hiT = 10**3.25, 10**4.5
loN, hiN = 10**-6.25, 10**-2.25
numbins = 200

dataloc = '/home/jacob/research/velas/vela2b/vela27/a0.490/'
fname = '{0:s}vela2b-27_GZa0.490.h5'.format(dataloc)

df = pd.read_hdf(fname, 'data')

index = ( (df['temperature']>loT) & (df['temperature']<hiT) & 
            (df['density']>loN) & (df['density']<hiN) &
            (df['x']>0) & (df['z']>0) & (np.abs(df['y'])<300) )

cloud = df[index]
cloudLoc = cloud[['x','y','z']]
cloudVel = cloud[['vx','vy','vz']]
print('Number of samples = {0:d}'.format(len(cloudLoc)))

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


# Determine the speed of each cell
cloudVel['speed'] = np.sqrt(cloud['vx']**2 + cloud['vy']**2 + cloud['vz']**2 )
loSpeed = cloudVel['speed'].min()
hiSpeed = cloudVel['speed'].max()


adotb = cloudVel['vx']*u0 + cloudVel['vy']*u1 + cloudVel['vz']*u2
bdotb = u0**2 + u1**2 + u2**2
factor = adotb/bdotb


cloudVel['along0'] = factor*u0
cloudVel['along1'] = factor*u1
cloudVel['along2'] = factor*u2
cloudVel['along'] = np.sqrt(cloudVel['along0']**2 + 
                    cloudVel['along1']**2 + cloudVel['along2']**2)

cloudVel['perp0'] = cloudVel['vx'] - cloudVel['along0']
cloudVel['perp1'] = cloudVel['vy'] - cloudVel['along1']
cloudVel['perp2'] = cloudVel['vz'] - cloudVel['along2']
cloudVel['perp'] = np.sqrt(cloudVel['perp0']**2 + 
                    cloudVel['perp1']**2 + cloudVel['perp2']**2)


# Plot historgram of dist
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
n, bins, patches = ax1.hist(cloudVel['speed'], bins=numbins, 
                            range=[loSpeed,hiSpeed],
                            histtype='step')
#n, bins, patches = ax2.hist(cloudVel['speed'], bins=numbins, 
#                            range=[loSpeed,hiSpeed],
#                            histtype='step',normed=True,cumulative=True)
n, bins, patches = ax2.hist(cloudVel['along'], bins=numbins, 
                            histtype='step',label='V Par')
n, bins, patches = ax2.hist(cloudVel['perp'], bins=numbins, 
                            histtype='step', label='V Perp')

ax1.set_xlabel('Speed')
ax2.set_xlabel('Speed')
ax2.legend()
s = 'vela2b-27_inflow_velocity_coherence.png'
fig.savefig(s, bbox_inches='tight', dpi=300)

sys.exit()

# Fit a rayleigh curve to the data
y = cloudVel['speed']

# Perform fits and 
x = np.linspace(y.min(), y.max(), 10000)


print('\nRayleigh Distribution')
param = st.rayleigh.fit(y)

y = st.rayleigh.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])
ax1.plot(x, y, label='Rayleigh')
results = st.ks_2samp(locM['dist'],y)
print('KS-Test: ',results)
#results = st.anderson_ksamp([locM['dist'],y])
print('Anderson: ',results)
y = st.rayleigh.cdf(x, *param[:-2], loc=param[-2], scale=param[-1])
ax2.plot(x, y, label='Rayleigh')
rayParam = param


#pd.DataFrame.hist(locM,column='dist',ax=ax,bins=numbins,color='r')
rayVar = st.rayleigh.var(loc=rayParam[-2],scale=rayParam[-1])
weiVar = st.weibull_max.var(c=weiParam[0], loc=weiParam[1],scale=weiParam[2])
distVar = locM['dist'].var()

print('Rayleigh Variance:     {0:.3f}'.format(rayVar))
print('Weibull Variance:      {0:.3f}'.format(weiVar))
print('Distribution Variance: {0:.3f}'.format(distVar))

# Compute standard deviation
rayStd = st.rayleigh.std(loc=rayParam[-2],scale=rayParam[-1])
weiStd = st.weibull_max.std(c=weiParam[0], loc=weiParam[1],scale=weiParam[2])
distStd = locM['dist'].std()

print('Rayleigh Std Dev:     {0:.3f}'.format(rayStd))
print('Weibull Std Dev:      {0:.3f}'.format(weiStd))
print('Distribution Std Dev: {0:.3f}'.format(distStd))


















