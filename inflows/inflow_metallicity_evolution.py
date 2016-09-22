
'''
Plots evolution of the mean mettalicity in each of the
narrow density bins
'''

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


loc = '/home/jacob/research/velas/vela2b/vela27/'
fname = loc+'a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'
rotname = loc+'a{0:.3f}/rotmat_a{0:.3f}.txt'



loT, hiT = 3.5, 4.5
cloudLimits = 'cloudLimits.csv'
limits = pd.read_csv(cloudLimits)

expns0 = range(20,50)
expns = [i/100. for i in expns0]

outfile = 'vela2b-27_inflow_metallicity_evolution.h5'
header = ['a','redshift','loN','hiN','loT','hiT','minSNII','maxSNII','meanSNII','medianSNII','sigdevSNII']
results = np.zeros(len(header))


for a in expns:

    print(a)
    nRange = limits[limits.expn==a]
    loN = nRange['loN'].values[0]
    hiN = nRange['hiN'].values[0]

    df = pd.read_hdf(fname.format(a), 'data')
    
    index = ( (df['temperature']>10**loT) & (df['temperature']<10**hiT) &
                (df['density']>10**loN) & (df['density']<10**hiN) &
                (df['x']<0) & (df['z']>0) & (np.abs(df['y'])<300) )

    cloud = df[index]

    res = np.zeros(len(header))

    res[0] = a
    res[1] = 1./a - 1
    res[2] = loN
    res[3] = hiN
    res[4] = loT
    res[5] = hiT
    res[6] = cloud['SNII'].min()
    res[7] = cloud['SNII'].max()
    res[8] = cloud['SNII'].mean()
    res[9] = cloud['SNII'].median()
    res[10] = cloud['SNII'].std()
    
    results = np.vstack((results,res))

results = np.delete(results,(0),axis=0)
df = pd.DataFrame(results,columns=header)
df.to_hdf(outfile,'data',mode='w')

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(results[:,0],results[:,6],label='Min')
ax.plot(results[:,0],results[:,8],label='Mean')
ax.plot(results[:,0],results[:,9],label='Median')
ax.plot(results[:,0],results[:,7],label='Max')
ax.set_yscale('log')
ax.legend(ncol=4,loc=9)
ax.set_xlabel('a')
ax.set_ylabel('SNII')
s = 'inflow_Z_evolution.png'
fig.savefig(s,bbox_inches='tight',dpi=300)












