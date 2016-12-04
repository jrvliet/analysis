
'''
Plots the phase diagram of the galaxy box. Produces three versions with seperate
statistics. One is colored by the mean SNII mass fraction in the bin, one is
colored by the minimum SNII mass fraction, one colored by the 25% quartile of
SNII mass fraction.
'''


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st


def mfToZ(mf):

    X_sun = 0.70683
    Y_sun = 0.27431
    Z_sun = 0.0188 
    r = 0.3347
    #X_cell = (1 - mf) / (1 + r)
    #Z = 
    return (mf / ((1 - mf) / (1 + r))) / (Z_sun / X_sun)


def quartile(a):
    return np.percentile(a,25)


def mkHist(ax,n,t,z,stat,cbarLabel):
    
    numbins = 200
    binrange = [[-10,2],[2,8]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(n,t,z,
                                statistic=stat,range=binrange,
                                bins=numbins)

    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel('Density [cm$^{-3}$]')
    ax.set_ylabel('Temperature [K]')
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(cbarLabel,rotation=270,fontsize=12)

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'

expns = np.arange(0.200,0.550,0.01)


for a in expns:

    print(a)
    redshift = 1./a - 1
    fname = dataloc+filename.format(a)

    df = pd.read_hdf(fname,'data')

    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
    cbarLabel = 'SNII Mass Fraction'
    mkHist(ax1,np.log10(df['density']),np.log10(df['temperature']),
            df['SNII'],'mean',cbarLabel)
    mkHist(ax2,np.log10(df['density']),np.log10(df['temperature']),
            df['SNII'],np.min,cbarLabel)
    mkHist(ax3,np.log10(df['density']),np.log10(df['temperature']),
            df['SNII'],np.max,cbarLabel)

    ax1.set_title('Mean')
    ax2.set_title('Min')
    ax3.set_title('Max')
    #ax3.set_title('25\% Quartile')

    #fig.suptitle('a={0:.3f}, z={1:.3f}'.format(a,redshift))
    fig.tight_layout()
    
    s = 'phase_SNII_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)

    df['metal'] = mfToZ(df['SNII']+df['SNIa'])
    
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    cbarLabel = 'log(Z)'
    mkHist(ax,np.log10(df['density']),np.log10(df['temperature']),
           df['metal'],np.min,cbarLabel) 

    s = 'initphase_min_a{0:.3f}.pdf'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


