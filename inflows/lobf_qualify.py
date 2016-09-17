
'''
Determines the quality of the line of best fit (LOBF) used
for the ''ideal'' inflow in vela2b-27

Methods:
1) Determine the impact parameter of the LOBF -> plot evolution
2) Plot components of LOBF -> plot evolution
3) Determine the LOBF for both full density range and narrow range
'''

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg as sl
import scipy.stats as st

pd.options.mode.chained_assignment = None


def fitLine(loc):
    
    locM = loc - loc.mean()
    covar = locM.cov()
    (u,s,v) = sl.svd(covar)
    u0 = u[0,0]
    u1 = u[1,0]
    u2 = u[2,0]
    
    return u0,u1,u2

def impact(u0,u1,u2,mean):

    x = mean['x']
    y = mean['y']
    z = mean['z']
    
    top2a = (u2*y - u1*z)**2
    top2b = (u0*z - u2*x)**2
    top2c = (u1*x - u0*y)**2
    top = np.sqrt( top2a + top2b + top2c) 
    
    bot = np.sqrt(u0**2 + u1**2 + u2**2)

    d = top/bot

    return d


# File paths
dataloc = '/home/jacob/research/velas/vela2b/vela27/'
basename = '{0:s}a{1:.3f}/vela2b-27_GZa{1:.3f}.h5'

# Full ranges
loT, hiT = 10**3.5, 10**4.5

# Snapshot expansion parameters
expns = np.arange(0.200,0.500,0.01)
expns0 = range(20,50)
expns = [i/100. for i in expns0]

# Read in cloud narrow density limits
limits = pd.read_csv('cloudLimits.csv')

# Output files
outfile = 'vela2b-27_cloud_LOBF.h5'
header = ['expn','fullb','fullu0','fullu1','fullu2',
            'narrowb','narrowu0','narrowu1','narrowu2']
fit = np.zeros(len(header))


fb, fu0, fu1, fu2 = [], [], [], []
nb, nu0, nu1, nu2 = [], [], [], []

for a in expns:
    print(a)
    redshift = 1./a - 1

    # Read in the gasbox
    fname = basename.format(dataloc,a)
    df = pd.read_hdf(fname, 'data')

    fitparams = np.zeros(len(header))
    fitparams[0] = a

    #############################################################3

    # Full cut
    loN, hiN = -6.25, -2.25
    index = ( (df['temperature']>loT) & (df['temperature']<hiT) & 
                (df['density']>10**loN) & (df['density']<10**hiN) &
                (df['x']>0) & (df['z']>0) & (np.abs(df['y'])<300) )
    cloud = df[index]
    cloudLoc = cloud[['x','y','z']]
    
    fullu0, fullu1, fullu2 = fitLine(cloudLoc)
    fullmean = cloudLoc.mean() 
    fullb = impact(fullu0,fullu1,fullu2,fullmean)

    fitparams[1] = fullb
    fitparams[2] = fullu0
    fitparams[3] = fullu1
    fitparams[4] = fullu2
    fb.append(fullb)
    fu0.append(fullu0)
    fu1.append(fullu1)
    fu2.append(fullu2)


    #############################################################3

    # Narrow cut
    nRange = limits[limits.expn==a]
    loN = nRange['loN'].values[0]
    hiN = nRange['hiN'].values[0]
    index = ( (df['temperature']>loT) & (df['temperature']<hiT) & 
                (df['density']>10**loN) & (df['density']<10**hiN) &
                (df['x']>0) & (df['z']>0) & (np.abs(df['y'])<300) )
    cloud = df[index]
    cloudLoc = cloud[['x','y','z']]
    
    narrowu0, narrowu1, narrowu2 = fitLine(cloudLoc)
    narrowmean = cloudLoc.mean() 
    narrowb = impact(narrowu0,narrowu1,narrowu2,narrowmean)
    
    fitparams[5] = narrowb
    fitparams[6] = narrowu0
    fitparams[7] = narrowu1
    fitparams[8] = narrowu2
    nb.append(narrowb)
    nu0.append(narrowu0)
    nu1.append(narrowu1)
    nu2.append(narrowu2)

    fit = np.vstack((fit,fitparams))

    #############################################################3

    # Plot the lines
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(1,1,1,projection='3d')

    t = np.linspace(-100,100,1000)
    x, y, z = [],[],[]
    for i in range(len(t)):
        x.append(fullmean['x']+t[i]*fullu0)
        y.append(fullmean['y']+t[i]*fullu1)
        z.append(fullmean['z']+t[i]*fullu2)
    ax.plot(x,y,z,'b')

    x, y, z = [],[],[]
    for i in range(len(t)):
        x.append(narrowmean['x']+t[i]*narrowu0)
        y.append(narrowmean['y']+t[i]*narrowu1)
        z.append(narrowmean['z']+t[i]*narrowu2)
    ax.plot(x,y,z,'r')

    ax.view_init(elev=30, azim=0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('a = {0:.3f}, z = {0:.3f}'.format(a,redshift))

    s = 'lobf_3d_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)


# Plot
fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))
axes = [ax1,ax2,ax3,ax4]

ax1.plot(expns,fb,label='Full')
ax1.plot(expns,nb,label='Narrow')
ax1.set_ylabel('Impact Paramter')

ax2.plot(expns,fu0,label='Full')
ax2.plot(expns,nu0,label='Narrow')
ax2.set_ylabel('u0')

ax3.plot(expns,fu1,label='Full')
ax3.plot(expns,nu1,label='Narrow')
ax3.set_ylabel('u1')

ax4.plot(expns,fu2,label='Full')
ax4.plot(expns,nu2,label='Narrow')
ax4.set_ylabel('u2')

for ax in axes:
    ax.set_xlabel('Expn')
    ax.legend()

fig.tight_layout()
s = 'lobf_qualify.png'
fig.savefig(s,bbox_inches='tight',dpi=300)


# Write the fit paramters to file
fit = np.delete(fit,(0),axis=0)
df = pd.DataFrame(fit,columns=header)
df.to_hdf(outfile,'data',mode='w')
























