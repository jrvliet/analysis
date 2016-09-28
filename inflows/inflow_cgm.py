
'''
Plots the spatial distribution of gas that fits the inflow
phase cut that lies within 1 Rvir of the galaxy
'''


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def mkPlot(x1,y1,z1,x2,y2,z2,a):

    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,3))
    ax1.scatter(x1,y1,marker='o',color='green',alpha=0.01)
    ax2.scatter(x1,z1,marker='o',color='green',alpha=0.01)
    ax3.scatter(y1,z1,marker='o',color='green',alpha=0.01)
    ax1.scatter(x2,y2,marker='o',color='plum',alpha=0.01)
    ax2.scatter(x2,z2,marker='o',color='plum',alpha=0.01)
    ax3.scatter(y2,z2,marker='o',color='plum',alpha=0.01)
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')

    ax2.set_xlabel('x')
    ax2.set_ylabel('z')

    ax3.set_xlabel('y')
    ax3.set_ylabel('z')

    ax1.set_xlim([-1,1])
    ax1.set_ylim([-1,1])
    ax2.set_xlim([-1,1])
    ax2.set_ylim([-1,1])
    ax3.set_xlim([-1,1])
    ax3.set_ylim([-1,1])

    redshift = 1./a - 1.
    ax2.set_title('a = {0:.3f}, z = {1:.3f}'.format(a,redshift))

    fig.tight_layout()
    savename = 'inflow_cgm_a{0:.3f}.png'.format(a)
    fig.savefig(savename,bbox_inches='tight',dpi=300)


baseloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotmat = 'a{0:.3f}/rotmat_a{0:.3f}.txt'
limitsFile = 'cloudLimits.csv'

limits = pd.read_csv(limitsFile)
loT, hiT = 10**3.5, 10**4.5

expns0 = range(20,50)
expns = [i/100. for i in expns0]

for a in expns:

    print(a)

    with open(baseloc+rotmat.format(a), 'r') as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = baseloc+filename.format(a)
    df = pd.read_hdf(fname, 'data')


    nRange = limits[limits.expn==a]
    loN = 10**(nRange['loN'].values[0])
    hiN = 10**(nRange['hiN'].values[0])
    
    index = ( (df['temperature']>loT) & (df['temperature']<hiT) &
                (df['density']>loN) & (df['density']<hiN) & df['r']<=rvir)

    
    index2 = index & (df['x']<0) & (df['z']>0)
    index1 = index & ~index2

    inflow = df[index1]
    cloud = df[index2]

    x1 = inflow['x']/rvir
    y1 = inflow['y']/rvir
    z1 = inflow['z']/rvir

    x2 = cloud['x']/rvir
    y2 = cloud['y']/rvir
    z2 = cloud['z']/rvir

    mkPlot(x1,y1,z1,x2,y2,z2,a)



    
                







































