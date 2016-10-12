


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


boxloc = '/mnt/cluster/cgm/vela2b/vela27/'
boxname = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotmat = 'a{0:.3f}/rotmat_a{0:.3f}.txt'
denseloc = '/lustre/projects/p089_swin/jvander/analysis/inflows/'
densename = 'projected_distance_distribution.h5'

dense = pd.read_hdf(denseloc+densename, 'data')
dense['std'] = np.sqrt(dense['speedStd']**2 + dense['stdDev']**2)

d = dense[dense['numcells']>1e4]
groups = d.groupby('a')

loT,hiT = 3.5,4.5

for a, group in groups:

    # Get the denstiy limits
    minIndex = group['std'].argmin()
    loN = group['loN'].ix[minIndex]
    hiN = group['hiN'].ix[minIndex]

    # Get Rvir
    with open(boxloc+rotmat.format(a)) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    # Read in the box
    fname = boxloc + boxname.format(a)
    df = pd.read_hdf(fname,'data')

    tempInds = (df['temperature']<10**hiT) & (df['temperature']>10**loT)
    densInds = (df['density']<10**hiN) & (df['density']>10**loN)
    spacInds = (df['zRot']>0) & (df['theta']<80)

    cloud = df[tempInds & densInds & spacInds]

    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
    
    ax1.scatter(cloud['xRot']/rvir,cloud['yRot']/rvir,marker='o',color='green',alpha=0.01)
    ax2.scatter(cloud['xRot']/rvir,cloud['zRot']/rvir,marker='o',color='green',alpha=0.01)
    ax3.scatter(cloud['yRot']/rvir,cloud['zRot']/rvir,marker='o',color='green',alpha=0.01)

    ax1.set_xlabel('x')
    ax2.set_xlabel('x')
    ax3.set_xlabel('y')
    ax1.set_ylabel('y')
    ax2.set_ylabel('z')
    ax3.set_ylabel('z')

    ax1.set_xlim([-3,3])
    ax2.set_xlim([-3,3])
    ax3.set_xlim([-3,3])
    ax1.set_ylim([-3,3])
    ax2.set_ylim([0,6])
    ax3.set_ylim([0,6])

    ax2.set_title('a={0:.3f}, z={1:.3f}, nH={2:.1f}-{3:.1f}'.format(
                    a,1./a-1,loN,hiN))

    fig.tight_layout()
    s = 'density_cut_spatialLoc_a{0:.3f}.png'.format(a)
    fig.savefig(a,bbox_inches='tight',dpi=300)
    plt.close(fig)











