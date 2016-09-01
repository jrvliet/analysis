
'''
Plots the radial velocity of cells in the inflows
of vela2b-27 as a funciton of galactocentric distance
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fname = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'


loT, hiT = 10**3, 10**4.5
loN, hiN = 10**-6, 10**-2.25


expns = np.arange(0.200,0.510,0.01)

for a in expns:

    z = 1./a - 1.0

    fn = fname.format(a)

    df = pd.read_hdf(fn,'data')

    index = ( (df['temperature']<hiT) & (df['temperature']>loT) &
                (df['density']<hiN) & (df['density']>loN) & 
                (df['x']<0) & (df['z']>0) & (np.abs(df['y'])<300) )

    cloud = df[index]

    cloud['r'] = np.sqrt(cloud['x']**2 + cloud['y']**2 + cloud['z']**2)
    cloud['vr'] = ( cloud['vx']*cloud['x'] + cloud['vy']*cloud['y'] + 
                    cloud['vz']*cloud['z'] ) / cloud['r']

    
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    ax.plot(cloud['r'], cloud['vr'], 'ok', alpha=0.01)
    ax.set_xlabel('Galactocentric Distance [kpc]')
    ax.set_ylabel('Radial Velocity [km/s]')
    ax.set_title('z = {0:.3f}'.format(z))

    s = 'vela2b-27_a{0:.3f}_cloud_radVel.png'.format(a)
    fig.savefig(s, bbox_inches='tight', dpi=250)

    
    










