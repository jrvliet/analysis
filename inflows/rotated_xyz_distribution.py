

'''
Plots xy, xz, yz distribution of points in rotated box
'''


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def mkPlot(x,y,z,a):

    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))

    ax1.scatter(x,y,marker='o',color='green',alpha=0.01)
    ax2.scatter(x,z,marker='o',color='green',alpha=0.01)
    ax3.scatter(y,z,marker='o',color='green',alpha=0.01)
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')

    ax2.set_xlabel('x')
    ax2.set_ylabel('z')

    ax3.set_xlabel('y')
    ax3.set_ylabel('z')

    ax1.set_xlim([-3,3])
    ax1.set_ylim([-3,3])
    ax2.set_xlim([-3,3])
    ax2.set_ylim([-3,3])
    ax3.set_xlim([-3,3])
    ax3.set_ylim([-3,3])

    redshift = 1./a - 1.
    ax2.set_title('a = {0:.3f}, z = {1:.3f}'.format(a,redshift))

    fig.tight_layout()
    savename = 'inflow_rotated_xyz_a{0:.3f}.png'.format(a)
    fig.savefig(savename,bbox_inches='tight',dpi=300)


# File parameters
baseloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotmat = 'a{0:.3f}/rotmat_a{0:.3f}.txt'
expns = np.arange(0.200,0.500,0.01)

for a in expns:
    print(a)
    
    with open(baseloc+rotmat.format(a)) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = baseloc+filename.format(a)
    df = pd.read_hdf(fname, 'data')
    
    x = df['xRot']/rvir
    y = df['yRot']/rvir
    z = df['zRot']/rvir

    mkPlot(x,y,z,a)
