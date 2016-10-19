
from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gc
import statsmodels as sm
import sys
import scipy.stats as st
import matplotlib.ticker as ticker


def fmt(x,pos):
    return '{0:.0f}'.format(x)

def mkHist(ax,x,y,xlabel,ylabel):

    numbins = 50
    stat = 'count'

    if 'z' in ylabel:
        binrange = [[-3,3],[0,6]]
    else:
        binrange = [[-3,3],[-3,3]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,y,
                                statistic=stat, bins=numbins,
                                range=binrange)

    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    
    mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Count')



mH = 1.6737e-24
pc2cm = 3.086e18
mSun = 1.989e33

vmin, vmax = -8, -1

galNums = [27]
expns = np.arange(0.200,0.500,0.01) 

# Temp range loop
minTemps = [3.5,4.5,5.5]
maxTemps = [4.5,5.5,6.5]
tempLabels = ['cool','warm','hot']

simples = [(-5.5,-5.0),(-5.0,-4.5),(-4.5,-4.0),
            (-4.0,-3.5),(-3.5,-3.0),(-3.0,-2.5)]


for a in expns:
    print('\ta = {0:.3f}'.format(a))

#    fig,((ax1,ax2,ax3,ax4,ax5,ax6),
#         (ax7,ax8,ax9,ax10,ax11,ax12),
#         (ax13,ax14,ax15,ax16,ax17,ax18)) = plt.subplots(3,6,figsize=(18,9))


    for i in range(len(minTemps)):
        fig,axes =  plt.subplots(3,6,figsize=(18,9))
        print('\t{0:s} gas'.format(tempLabels[i]))
        loT = 10**minTemps[i]
        hiT = 10**maxTemps[i]
        print(minTemps[i],maxTemps[i])

        dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/'.format(a)
        dataloc = '/home/jacob/research/velas/vela2b/vela27/a{0:.3f}/'.format(a)
        boxfile = '{0:s}vela2b-27_GZa{1:.3f}.rot.h5'.format(dataloc,a)


        rotmatFile = '{0:s}rotmat_a{1:.3f}.txt'.format(dataloc,a)
        with open(rotmatFile) as f:
            f.readline()
            rvir = float(f.readline().split()[3])
            
        df = pd.read_hdf(boxfile, 'data')

        baseInds = ( (df['temperature']<=hiT) & (df['temperature']>=loT) & 
                     (df['zRot']>0) & (df['theta']<80) )

        for j in range(len(simples)):
            
            loN, hiN = simples[j][0], simples[j][1]     
            denseInds = (df['density']>10**loN) & (df['density']<10**hiN)

            cloud = df[baseInds & denseInds]
            print(loN,hiN,loT,hiT,len(cloud))

            mkHist(axes[0][j],cloud['xRot']/rvir,cloud['yRot']/rvir,'x','y')
            mkHist(axes[1][j],cloud['xRot']/rvir,cloud['zRot']/rvir,'x','z')
            mkHist(axes[2][j],cloud['yRot']/rvir,cloud['zRot']/rvir,'y','z')

            axes[0][j].set_title('{0:.1f} $<$ nH $<$ {1:.1f}'.format(loN,hiN))

        

        fig.tight_layout()
        s = 'density_cut_locations_a{0:.3f}_{1:s}.png'.format(a,tempLabels[i])
        fig.savefig(s,bbox_inches='tight',dpi=300)
        plt.close(fig)






