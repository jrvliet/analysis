
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
import scipy.stats as st


def mkPlot(x1,y1,z1,x2,y2,z2,a):


    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(15,10))

    mkHist(ax1,x1,y1,'x','y')
    mkHist(ax2,x1,z1,'x','z')
    mkHist(ax3,y1,z1,'y','z')

    mkHist(ax4,x2,y2,'x','y')
    mkHist(ax5,x2,z2,'x','z')
    mkHist(ax6,y2,z2,'y','z')


    for ax in [ax1,ax2,ax3]:
        ax.set_title('Inflow')
    for ax in [ax4,ax5,ax6]:
        ax.set_title('Cloud')

    redshift = 1./a - 1.
    ax2.set_title('a = {0:.3f}, z = {1:.3f}'.format(a,redshift))

    fig.tight_layout()
    savename = 'inflow_cgm_a{0:.3f}.png'.format(a)
    fig.savefig(savename,bbox_inches='tight',dpi=300)


def mkHist(ax,x,y,xlab,ylab):

    stat = 'count'
    numbins = 50
    binrange = [[-1,1],[-1,1]]

    h,xedges,yedges,binnumber = st.binned_statistic_2d(x1,y1,z1,
                                statistic=stat, bins=numbins,
                                range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    mesh = ax.pcolormesh(xedges, yedges, h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Counts',rotation=270,fontsize=12)
    

baseloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
baseloc = '/home/jacob/research/velas//vela2b/vela27/'
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

    print('Inflow size = {0:d}'.format(len(inflow['x'])))
    print('Cloud size  = {0:d}'.format(len(cloud['x'])))

    x1 = inflow['x']/rvir
    y1 = inflow['y']/rvir
    z1 = inflow['z']/rvir

    x2 = cloud['x']/rvir
    y2 = cloud['y']/rvir
    z2 = cloud['z']/rvir

    mkPlot(x1,y1,z1,x2,y2,z2,a)




    
                







































