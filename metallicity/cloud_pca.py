
'''
Selects gas from the simulations that fit the cloud phase
Computes the PCA analysis and writes the results to file
'''

import pandas as pd
import scipy.linalg as sl
import numpy as np
import matplotlib.pyplot as plt

dataloc = '/home/jacob/research/velas/vela2b/vela27/a0.490/'
boxfile = '{0:s}vela2b-27_GZa0.490.h5'.format(dataloc)

# Read in data
d = pd.read_hdf(boxfile, 'data')

# Output file:
f = open('cloud_size_pca.dat', 'w')
f.write('hiN\tNum\ts1\tv1\t\t\ts2\tv2\t\t\ts3\tv3\n')
string = ('{0:.2f}\t{13:d}\t{1:.1e}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t'
    '{5:.1e}\t{6:.4f}\t{7:.4f}\t{8:.4f}\t'
    '{9:.1e}\t{10:.4f}\t{11:.4f}\t{12:.4f}\n')



# Define cloud parameters
loT, hiT = 10**4, 10**4.5
loN = 10**-5.1

ns = np.linspace(-5,2,31)

eVal = []
number = []

for n in ns:
    print('hiN = 10**{0:f}'.format(n))
    hiN = 10**n

    # Select out gas that fits this paramters
    cloudInds = ( (d['temperature']<hiT) & (d['temperature']>loT) & 
                  (d['density']<hiN) & (d['density']>loN) &
                  (d['x']<0) & (d['z']>0) & (np.abs(d['y'])<300))
    cloud = d[cloudInds]

    # Only need the coordinates of the cells
    loc = cloud[['x','y','z']]

    # Normalize the cloud locations
    locM = loc - loc.mean()

    # Compute the covariance matrix
    covar = locM.cov()

    # Compute the single value decomposition of the 
    # covariance matrix. 
    # s = eignenvalues
    # u = eigenvectors
    (u,s,v) = sl.svd(covar)

    f.write(string.format(n, s[0], u[0,0], u[1,0], u[2,0],
                             s[1], u[0,1], u[1,1], u[2,1],
                             s[2], u[0,2], u[1,2], u[2,2], len(cloud)))

    eVal.append(s[0])
    number.append(len(cloud))
f.close()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(ns,eVal,'xk',label='Eigenvals')
ax2.plot(ns,number,'ob',label='Number')
ax1.set_ylabel('Largest Eigenvalue')
ax2.set_ylabel('Number of Cells')
plt.legend(frameon=False,loc='lower right')
plt.savefig('cloud_size_evals.png',bbox_inches='tight',dpi=300)




