
'''
Plots the radial velocity of cells in the inflows
of vela2b-27 as a funciton of galactocentric distance
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


fname = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'
rname = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/rotmat_a{0:.3f}.txt'


loT, hiT = 10**3, 10**4.5
loN, hiN = 10**-6, 10**-2.25

expns = np.arange(0.200,0.510,0.01)

minv, maxv = [],[]

for a in expns:

    z = 1./a - 1.0

    with open(rname.format(a),'r') as f:
        f.readline()
        rvir = float( f.readline().split()[3] ) 

    fn = fname.format(a)

    try:
        df = pd.read_hdf(fn,'data')
    except IOError:
        continue

    index = ( (df['temperature']<hiT) & (df['temperature']>loT) &
                (df['density']<hiN) & (df['density']>loN) & 
                (df['x']<0) & (df['z']>0) & (np.abs(df['y'])<300) )

    cloud = df[index]

    cloud['r'] = np.sqrt(cloud['x']**2 + cloud['y']**2 + cloud['z']**2)
    cloud['vr'] = ( cloud['vx']*cloud['x'] + cloud['vy']*cloud['y'] + 
                    cloud['vz']*cloud['z'] ) / cloud['r']

    posInds = (cloud['vr']>0)
    negInds = (cloud['vr']<0)
    poscloud = cloud[posInds]
    negcloud = cloud[negInds]
    minv.append(cloud['vr'].min())
    maxv.append(cloud['vr'].max())
    print(cloud['r'].max()/rvir)
    
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    
    ax.plot(poscloud['r']/rvir, poscloud['vr'], 'ob', alpha=0.01)
    ax.plot(negcloud['r']/rvir, negcloud['vr'], 'or', alpha=0.01)
    ax.set_xlim([0,5.0])
    ax.set_ylim([-500,500])
    ax.set_xlabel('Galactocentric Distance [kpc]')
    ax.set_ylabel('Radial Velocity [km/s]')
    ax.set_title('z = {0:.3f}'.format(z))

    s = 'vela2b-27_a{0:.3f}_cloud_radVel.png'.format(a)
    fig.savefig(s, bbox_inches='tight', dpi=250)

    
    
print( 'Min Rad Vel = {0:.2f}'.format(min(minv)))
print( 'Max Rad Vel = {0:.2f}'.format(max(maxv)))










