
'''
Determines the coherence of the inflows
Focus on the poster child inflow in vela2b-27
Outputs variance of data and scale of fit 
for all times
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
plotting = 0

loT, hiT = 10**3.25, 10**4.5
loN, hiN = 10**-6.25, 10**-2.25
numbins = 200

dataloc = '/home/jacob/research/velas/vela2b/vela27/a{0:.3f}/'
dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'

expns = np.arange(0.200,0.500,0.01)

outfilename = 'vela2b-27_cloud1_rayleightFit.dat'
fout = open(outfilename, 'w')
header = 'Expn\tNumber\tScale\tLoc\tRayVary\tDataVar\tKS\n'
fout.write(header)
form = '{0:.3f}\t{1:d}\t{2:.3f}\t{3:.3f}\t{4:.3f}\t{5:.3f}\t{6:.3f}\n'

for a in expns:
    print(a)

    fname = '{0:s}a{1:.3f}/vela2b-27_GZa{1:.3f}.h5'.format(dataloc,a)

    df = pd.read_hdf(fname, 'data')

    index = ( (df['temperature']>loT) & (df['temperature']<hiT) & 
                (df['density']>loN) & (df['density']<hiN) &
                (df['x']>0) & (df['z']>0) & (np.abs(df['y'])<300) )

    cloud = df[index]
    cloudLoc = cloud[['x','y','z']]
    numCells = len(cloudLoc)

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
    locM['dist'] = np.sqrt(locM['a0']**2 + locM['a1']**2 + locM['a2']**2)

    out = pd.cut(locM['dist'], bins=numbins)
    counts = pd.value_counts(out, sort=False)
    counts.to_csv('counts.out',sep=' ',index=True,header=False)
    normedCounts = counts/sum(counts)

    # Fit a rayleigh curve to the data
    y = locM['dist']

    # Perform fits and 
    x = np.linspace(y.min(), y.max(), 10000)

    param = st.rayleigh.fit(locM['dist'])
    y = st.rayleigh.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])
    results = st.ks_2samp(locM['dist'],y)
    ks = results[0]

    scale = param[-1]
    loc = param[-2]

    # Compute variance
    rayVar = st.rayleigh.var(loc=loc,scale=scale)
    distVar = locM['dist'].var()

    fout.write(form.format(a,numCells,scale,loc,rayVar,distVar,ks))


    # Plot
    fig, ax = plt.subplots(1,1,figsize=(5,5))
    ax.hist(locM['dist'], bins=numbins, range=[locM['dist'].min(),locM['dist'].max()],
            histtype='step', normed=True, label='Data')
    ax.plot(x,y,label='Rayleigh Fit')
    ax.set_xlabel('Distance from Line')
    redshift = 1./a - 1.
    ax.set_title('z = {0:.3f}'.format(redshift))

    s = './spatialCoherence/vela2b-27_a{0:.3f}_spatialCoherence.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)





fout.close()
















