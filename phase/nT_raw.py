#!/usr/bin/python
# Filename: nT_raw.py
#
# Author: Jacob Vander Vliet
# Version: 1
# Date: 4/12/2013
# Last Modified: 4/12/13
#
# Description:
#   Creates the phase diagram of the raw snapshots
#   
#
# Instructions:
#   python nT_bulk.py <normalize?> <plot together?>
#

import numpy as np
import matplotlib.pyplot as plt
from os import system
import scipy.interpolate
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib as mpl
import sys

params = {'xtick.labelsize':11,
          'ytick.labelsize':11,
          'font.size':12}
plt.rcParams.update(params)


galID_list = ['D9o2', 'D9q', 'D9m4a']
feedback_list = ['dwSN', 'dwALL_1', 'dwALL_8']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]


file_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/'
numbins = 50
normalize = 0
histos = []
xed = []
yed = []

# Create the subplots
fig,((p11,p12,p13)) = plt.subplots(1,3,figsize=(10.2, 3.5))
plot_list = [p11, p12, p13]

for i in range(0,len(galID_list)):

    galID = galID_list[i]
    expn = expn_list[i]

    n = []
    t = []

    for a in expn:

        # Open the data
        abs_file = galID+'_outputs/a'+a+'/'+galID+'_GZa'+a+'.txt'
        lognH, logT = np.loadtxt(file_loc+abs_file, skiprows=2, usecols=(7, 8), unpack=True)
        print galID, a

        for j in range(0,len(lognH)):

            n.append(np.log10(lognH[j]))
            t.append(np.log10(logT[j]))

    H, xedges, yedges = np.histogram2d( n, t, bins=numbins, normed=normalize)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    print 'Density: ',min(n),max(n)
    print 'Temperature: ',min(t), max(t)
    # Rotate and filp the histogram
    H = np.rot90(H)
    H = np.flipud(H)
    
    # Mask the bins where the count is zero
    Hmasked = np.ma.masked_where(H==0,H)
        
    # Take the log of the count
    Hmasked = np.log10(Hmasked)
    xed.append(xedges)
    yed.append(yedges)
    histos.append(Hmasked)



labels = ['(a)', '(b)', '(c)']
print ''
print len(histos)
print ''
for i in range(0,len(plot_list)):
    ax = plot_list[i]
    H = histos[i]
    xedges = xed[i]
    yedges = yed[i]
    mesh = ax.pcolormesh(xedges,yedges,H)
    ax.set_xlim([-8, 1])
    ax.set_ylim([2,8])
    ax.text(-1, 7, labels[i])

    ax.set_ylabel('$\log$ (T [K] )')
    ax.set_xlabel(' $\log (n_{H} $ [cm$^{-3}$] )')

    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.set_label('$\log$ (Counts)')
#    plt.savefig('dum_'+str(i)+'.eps')
    
#dropx = [p12, p22, p32, p13, p23, p33, p11, p21, p31]
#plt.setp([a.get_xticklabels() for a in dropx],visible=False)

plt.tight_layout()
s = file_loc+'phase_bulk.eps'
plt.savefig(s, bbox_inches='tight')
