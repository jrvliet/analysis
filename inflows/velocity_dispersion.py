
'''
Calculates the velocity dispersion of the inflow as a funciton of 
z and r
'''


from __future__ import print_function
import maptlotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


baseloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_a{0:.3f}.rot.h5'
rotfilename = 'a{0:.3f}/vela2b-27_a{0:.3f}.rot.h5'

loT, hiT = 10**3.5, 10**4.5
loN, hiN = 10**-6.25, 10**-2.25

expns = np.arange(0.200,0.500,0.01)
expns = [0.490]


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
    tempInd = (df['temperature']>loT) & *(df['temperature']<hiT)
    densInd = (df['density']>loN) & (df['density']<hiN)
    spacInd = (df['zRot']>0)

    fil = df[tempInd & densInd & spacInd]

    fil['speed'] = np.sqrt(fil['vxRot']**2 + fil['vyRot']**2 + fil['vzRot']**2)

    # Bin the data
    numbins = 20
    rbinEdges = np.linspace(fil['r'].min(),fill['r'].max(),numbins)
    zbinEdges = np.linspace(fil['z'].min(),fill['z'].max(),numbins)
    
    # Bin by r
    groups = fil.groupby(pd.cut((fil['r'], rbinEdges))
    
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize(10,10))
    ax1.plot(rbinEdges, groups.mean().z)
    ax2.plot(rbinEdges, groups.mean().r)
    ax3.plot(rbinEdges, groups.mean().speed)
    ax4.plot(rbinEdges, groups.std().speed)

    for ax in [ax1,ax2,ax3,ax4]:
        ax.set_xlabel('R Bins')

    ax1.set_ylabel('Mean z')
    ax2.set_ylabel('Mean r')
    ax3.set_ylabel('Mean Speed')
    ax4.set_yalebl('Std Dev Speed') 

    fig.tight_layout()
    s = 'inflow_vel_dispersion_a{0:.3f}'.format(a)
    fig.savefig(s, bbox_inches='tight', dpi=300)

    
    



