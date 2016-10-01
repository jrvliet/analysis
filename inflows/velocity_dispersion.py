
'''
Calculates the velocity dispersion of the inflow as a funciton of 
z and r
'''


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


baseloc = '/home/jacob/research/velas/vela2b/vela27/'
baseloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotfilename = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

loT, hiT = 10**3.5, 10**4.5
loN, hiN = 10**-6.25, 10**-2.25

expns = [0.490]
expns = np.arange(0.200,0.500,0.01)


for a in expns:

    print(a)

    # Get Rvir
    rotfile = baseloc+rotfilename.format(a)
    with open(rotfile,'r') as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    # Read in box
    fname = baseloc+filename.format(a)
    df = pd.read_hdf(fname, 'data')

    # Select out filament
    tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
    densInd = (df['density']>loN) & (df['density']<hiN)
    spacInd = (df['zRot']>0)

    fil = df[tempInd & densInd & spacInd]

    fil['speed'] = np.sqrt(fil['vxRot']**2 + fil['vyRot']**2 + fil['vzRot']**2)

    # Bin the data
    numbins = 50
    fil['rmod'] = fil['r']/rvir
    fil['zmod'] = fil['z']/rvir
    rbinEdges = np.linspace(fil['rmod'].min(),fil['rmod'].max(),numbins)
    zbinEdges = np.linspace(fil['zmod'].min(),fil['zmod'].max(),numbins)

    # Bin by r
    groups = fil.groupby(pd.cut(fil['rmod'], rbinEdges))

    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
    
    ax1.plot(rbinEdges[:-1], groups.mean().speed)
    ax2.plot(rbinEdges[:-1], groups.std().speed)
    ax3.plot(rbinEdges[:-1], groups.size(),label='{0:.3f}'.format(a))

    for ax in [ax1,ax2,ax3]:
        ax.set_xlabel('R [Rvir]')

    ax1.set_ylabel('Mean Speed')
    ax2.set_ylabel('Std Dev Speed') 
    ax3.set_ylabel('Number of cells in group')

    fig.tight_layout()
    s = 'inflow_vel_dispersion_a{0:.3f}.png'.format(a)
    fig.savefig(s, bbox_inches='tight', dpi=300)

    
    



