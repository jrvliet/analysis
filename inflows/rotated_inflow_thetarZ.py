
'''
Plots a 2d histogram of the mean metals in the box
as a function of r, theta in the rotated gasbox.
The box is rotated so the inflow lies along the positive 
z-axis (theta=0)
'''

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st


def mkHist(ax,x,y,z,stat,numbins,binrange,xlab,ylab,clab,title):

    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,z,
                                    statistic=stat,bins=numbins,
                                    range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)

    mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(clab,rotation=270,fontsize=12)



# File parameters
baseloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotmat = 'a{0:.3f}/rotmat_a{0:.3f}.txt'
expns = np.arange(0.200,0.500,0.01)

# Plotting parameters
stat = 'mean'
numbins = 500
binrange = [[0,180],[0,6]]
xlab = 'Theta [deg]'
ylab = 'r [Rvir]'
clab = 'Mean SNII'
title = 'a = {0:.3f}, z = {1:.3f}'
savename = 'inflow_r_theta_SNII_a{0:.3f}.png'

for a in expns:
    print(a)
    z = 1./a - 1.

    # Get Rvir
    with open(baseloc+rotmat.format(a)) as f:
        f.readline()
        rvir = float(f.readline().split()[3])
    
    # Read in data
    fname = baseloc+filename.format(a)
    df = pd.read_hdf(fname, 'data')

    # Plot
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    mkHist(ax,df['theta'],df['r']/rvir,df['SNII'],stat,numbins,
            binrange,xlab,ylab,clab,title.format(a,z))
    fig.savefig(savename.format(a),bbox_inches='tight',dpi=300)

    
    
    
    




