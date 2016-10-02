
'''
Fits a 2d gaussian to the xy-distribuiton of
cell that fit the inflow phase as a funciton of z
'''

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



baseloc = '/home/jacob/research/vela2b/vela27/rotatedBox/'
baseloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotmatname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

# Inflow parameters
loT, hiT = 10**3.5, 10**4.5
loN, hiN = 10**-6.25, 10**-2.25

expns = ['0.490']
expns = np.arange(0.200,0.500,0.05)

fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))


for a in expns:

    print(a)

    # Get rvir
    with open(baseloc+rotmatname.format(a), 'r') as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    # Read box
    fname = baseloc+filename.format(a)
    df = pd.read_hdf(fname, 'data')
    
    tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
    densInd = (df['density']>loN) & (df['density']<hiN)
    spacInd = (df['zRot']>0)
    
    fil = df[tempInd & densInd & spacInd]
    
    # Bin the data
    numbins = 20
    zbinEdges = np.linspace(fil['zRot'].min(), fil['zRot'].max(), numbins)

    groups = fil.groupby(pd.cut( fil['zRot'], zbinEdges ))

    sx, sy, s, midz = [], [], [], []
    for key, group in groups:
        sigx = group['xRot'].std()
        sigy = group['yRot'].std()
        sigAll = np.sqrt( sigx**2 + sigy**2 )
        mid = group['zRot'].median()/rvir
    
        sx.append(sigx)
        sy.append(sigy)
        s.append(sigAll)
        midz.append(mid)
    
    ax1.plot(midz,sx,label=a)
    ax2.plot(midz,sy,label=a)
    ax3.plot(midz,s,label=a)

ax3.legend()
for ax in [ax1,ax2,ax3]:
    ax.set_xlabel('Z [Rvir]')
ax1.set_ylabel('Y Dispersion')
ax2.set_ylabel('X Dispersion')
ax3.set_ylabel('Total Dispersion')

fig.tight_layout()
s = 'spatial_dispersion_rotated.png'
fig.savefig(s,dpi=300,bbox_inches='tight')



