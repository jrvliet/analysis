
'''
Plots the radial velocity of cells in the inflows
of vela2b-27 as a funciton of galactocentric distance
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

base = '/mnt/cluster/abs/cgm/'
base = '/home/jacob/research/velas/'

fname = base+'vela2b/vela27/a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'
rname = base+'vela2b/vela27/a{0:.3f}/rotmat_a{0:.3f}.txt'
sname = base+'vela2b/vela27/subhalos/vela2b-27_cloud1_subhalos.h5'


subhalos = pd.read_hdf(sname,'data')


loT, hiT = 10**3, 10**4.5
loN, hiN = 10**-6, 10**-2.25

expns = np.arange(0.200,0.550,0.01)

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

    cgm = cloud[cloud['r']>0.1*rvir]

    cgm['vr'] = ( cgm['vx']*cgm['x'] + cgm['vy']*cgm['y'] + 
                    cgm['vz']*cgm['z'] ) / cgm['r']

    posInds = (cgm['vr']>0)
    negInds = (cgm['vr']<0)
    poscloud = cgm[posInds]
    negcloud = cgm[negInds]
    minv.append(cgm['vr'].min())
    maxv.append(cgm['vr'].max())
    
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    
    ax.plot(poscloud['r']/rvir, poscloud['vr'], 'o', color='cyan', alpha=0.01)
    ax.plot(negcloud['r']/rvir, negcloud['vr'], 'o', color='salmon', alpha=0.01)


    # Plot subhalos
    subs = subhalos[subhalos['a']==a]
    for i in range(len(subs)):
        d = subs['r'].iloc[i]
        vr = subs['vr'].iloc[i]
        subRvir = subs['rvir'].iloc[i]
        radius = subRvir/rvir
    
        print(a,z,d,vr,subRvir,radius)

        #circle = plt.Circle((d,vr), radius, color='k')
        #ax.add_artist(circle)
        ax.vlines(d,-500,500,colors='purple')
        lab = '{0:.2f}'.format(subs['MR'].iloc[i])
        ax.text(d,300,lab,fontsize=20)


    ax.set_xlim([0,5.0])
    ax.set_ylim([-500,500])
    ax.set_xlabel('Galactocentric Distance [kpc]')
    ax.set_ylabel('Radial Velocity [km/s]')
    ax.set_title('z = {0:.3f}'.format(z))

    s = 'vela2b-27_a{0:.3f}_cloud_radVel.png'.format(a)
    fig.savefig(s, bbox_inches='tight', dpi=200)


    
    
print( 'Min Rad Vel = {0:.2f}'.format(min(minv)))
print( 'Max Rad Vel = {0:.2f}'.format(max(maxv)))










