
'''
Plots the evolution of cloud 1 in vela2b-27
Points are colored by the density
'''

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


datafile = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'
rotmatfile = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/rotmat_a{0:.3f}.txt'

loT,hiT = 10**3.25, 10**4.5
loN,hiN = 10**-6.25,10**-2.25

xmin, xmax = [], []
ymin, ymax = [], []
zmin, zmax = [], []

expns = np.arange(0.200,0.500,0.010)

for a in expns:
    
    z = 1./a -1
    print('a = {0:.3f}, z = {1:.3f}'.format(a,z))

    with open(rotmatfile.format(a), 'r') as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    df = pd.read_hdf(datafile.format(a), 'data')

    index = ( (df['temperature']<hiT) & (df['temperature']>loT) &
                (df['density']<hiN) & (df['density']>loN) &
                (df['x']<0) & (df['z']>0) & (np.abs(df['y'])<300) )

    cloud = df[index]

    xmin.append( df['x'].min()/rvir )
    xmax.append( df['x'].max()/rvir )
    ymin.append( df['y'].min()/rvir )
    ymax.append( df['y'].max()/rvir )
    zmin.append( df['z'].min()/rvir )
    zmax.append( df['z'].max()/rvir )

    
    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(9,3))
    
    ax1.scatter(cloud['x']/rvir, cloud['y']/rvir, marker='o', c=cloud['SNII'].values, alpha=0.01)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_xlim([0,-3])
    ax1.set_ylim([-3,3])

    ax2.scatter(cloud['x']/rvir, cloud['z']/rvir, marker='o', c=cloud['SNII'].values, alpha=0.01)
    ax2.set_xlabel('x')
    ax2.set_ylabel('z')
    ax2.set_xlim([0,-3])
    ax2.set_ylim([0,3])
    
    ax3.scatter(cloud['y']/rvir, cloud['z']/rvir, marker='o', c=cloud['SNII'].values, alpha=0.01)
    ax3.set_xlabel('y')
    ax3.set_ylabel('z')
    ax3.set_xlim([-3,3])
    ax3.set_ylim([0,3])


    ax2.set_title('z = {0:.3f}'.format(z))

    fig.tight_layout()
    s = 'vela2b-27_a{0:.3f}_cloud1_densityColor.png'.format(a)
    fig.savefig(s.format(a), bbox_inches='tight', dpi=300)

    
        
print('X Limits = {0:.3f} - {1:.3f}'.format(min(xmin),max(xmax)))
print('Y Limits = {0:.3f} - {1:.3f}'.format(min(ymin),max(ymax)))
print('Z Limits = {0:.3f} - {1:.3f}'.format(min(zmin),max(zmax)))




