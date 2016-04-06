'''
Plots phase space of LOS position and LOS velocity
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.stats.mstats import gmean
import sys

def plotHist(hs, xeds, yeds, ions, filename):

    # Plot
    fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(1,4,figsize=(20.8, 5.2))
    axes = [ax1, ax2, ax3, ax4]

    for i in range(0,len(axes)):
        ion = ions[i]
        ax = axes[i]
        H = hs[i]
#        H = np.ma.masked_where(H==0,H)
#        H = np.log10(H)
        xedges = xeds[i]
        yedges = yeds[i]

        cm = mp.cm.get_cmap('viridis')
        mesh = ax.pcolormesh(xedges,yedges,H, cmap=cm)
        ax.set_ylabel('Dist Along LOS')
        ax.set_xlabel('LOS Velocity')
        ax.set_title(ion)
        if 'normed' in filename:
            ax.set_ylim([0,1])
        ax.set_xlim([-100, 100])
        ax.set_ylim((yedges[0],yedges[-1]))

        cbarLabel = '$\log$ (Geo Mean of $N_{ion}$)'
        cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel(cbarLabel, rotation=270, fontsize=12)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    fig.savefig(filename, bbox_inches='tight')



def make_histos(v, s, nIon, size, nbins, ion):

    '''
    Bins up the x and y data into nbins
    Value is each bin is the geometric mean of column density
    '''
    
    # Find column density
    column = np.zeros(len(size))
    for i, (n,l) in enumerate(zip(nIon, size)):
        # Take the cube root of the cell length and convert from kpc to cm
        length = l**(1./3.) * 3.086e21
        column[i] = (10**n) * length
    
    print min(column)
    print max(column)
    vmin = -100
    vmax = 100
    smin = -220
    smax = 220  

    print 'Making histogram'
    H, xed, yed = np.histogram2d(v,s,bins=nbins, range=[[vmin,vmax],[smin,smax]])
    h = np.zeros_like(H)
    print 'Histogram done\n'

    for i in range(0,H.shape[0]):
        for j in range(0,H.shape[1]):
            # #rows = shape[0]
            # #cols = shape[1]
            vmin = xed[j]
            vmax = xed[j+1]
            smin = yed[i]
            smax = yed[i+1]
            val = []
            for k, (vel,sloc) in enumerate(zip(v,s)):
                if vel>vmin and vel<vmax and sloc>smin and sloc<smax:
                    val.append(column[k])
            print np.log10(min(val))
            print np.log10(max(val))
            print np.log10(np.mean(val))
            print np.log10(gmean(val))
            print vmin, vmax
            print smin, smax
            h[i,j] = gmean(val)
            print h[i,j]
    h = np.log10(h)
    np.savetxt('{0:s}_velHist.out'.format(ion), h)
    return h,xed,yed


ions = ['HI', 'MgII', 'CIV', 'OVI']
numbins = 50

histos1, xed1, yed1 = [], [], []
histos2, xed2, yed2 = [], [], []
for ion in ions:
    
    print ion
    # Read in data    
    filename = '{0:s}_vlos.dat'.format(ion)
    l, s, v, n, t, Z, size, nIon = np.loadtxt(filename, skiprows=1, usecols=(0,1,2,3,4,5,6,7), unpack=True)
   
    mid = max(l) / 2.0
    for i in range(0,len(s)):
        s[i] = s[i]*l[i] - mid
    
    print min(s), max(s)
    # Make histogram
    H, xedges, yedges = make_histos(v, s, nIon, size, numbins, ion)

    # Make histogram for non-normed s
#    H, xedges, yedges = np.histogram2d(v,s, bins=numbins)
    H = np.rot90(H)
    H = np.flipud(H)
    xed2.append(xedges)
    yed2.append(yedges)
    histos2.append(H)

#plotHist(histos1, xed1, yed1, ions, 'vel2dHist_normed.png')
plotHist(histos2, xed2, yed2, ions, 'vel2dHist_gmeanColDense.pdf')









