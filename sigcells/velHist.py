'''
Plots phase space of LOS position and LOS velocity
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp

def plotHist(hs, xeds, yeds, ions, filename):

    # Plot
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(11.2, 10.8))
    axes = [ax1, ax2, ax3, ax4]

    for i in range(0,len(axes)):
        ion = ions[i]
        ax = axes[i]
        H = hs[i]
        H = np.ma.masked_where(H==0,H)
        H = np.log10(H)
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

        cbarLabel = '$\log$ (Counts)'
        cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel(cbarLabel, rotation=270, fontsize=12)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    fig.savefig(filename, bbox_inches='tight')



ions = ['HI', 'MgII', 'CIV', 'OVI']
numbins = 50

histos1, xed1, yed1 = [], [], []
histos2, xed2, yed2 = [], [], []
for ion in ions:
    
    # Read in data    
    filename = '{0:s}_vlos.dat'.format(ion)
    l, s, v = np.loadtxt(filename, skiprows=1, usecols=(0,1,2), unpack=True)
   
    print min(s)
    print max(s)    
 
    # Make histogram
    H, xedges, yedges = np.histogram2d(v,s, bins=numbins)
    H = np.rot90(H)
    H = np.flipud(H)
    xed1.append(xedges)
    yed1.append(yedges)
    histos1.append(H)

    # Make histogram for non-normed s
    for i in range(0,len(s)):
        s[i] = s[i]*l[i]
    H, xedges, yedges = np.histogram2d(v,s, bins=numbins)
    H = np.rot90(H)
    H = np.flipud(H)
    xed2.append(xedges)
    yed2.append(yedges)
    histos2.append(H)

plotHist(histos1, xed1, yed1, ions, 'vel2dHist_normed.png')
plotHist(histos2, xed2, yed2, ions, 'vel2dHist.png')

