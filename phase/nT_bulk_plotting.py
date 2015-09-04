#!/usr/bin/python
# Filename: nT_bulk.py
#
# Author: Jacob Vander Vliet
# Version: 2
# Date: 4/12/2013
# Last Modified: 4/12/13
#
# Description:
#   Creates the phase diagram of the gas for cells that
#   contribute to absorption
#
# Instructions:
#   python nT_bulk.py <normalize?>
#
# If normalize is high, the ALL1 and ALL8 are normalized to the SN counts

import numpy as np
import matplotlib.pyplot as plt
from os import system
import scipy.interpolate
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib as mpl
import sys


if len(sys.argv)==1:
    print 'Usage:\n\t python nT_bulk_plotting.py <normalize?>\n'
    sys.exit()

verbose = 0
normalize = int(sys.argv[1])
bulk = 1

galID_list = ['D9o2', 'D9q', 'D9m4a']
feedback_list = ['dwSN', 'dwALL_1', 'dwALL_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']


file_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/abscells/'
numbins = 50
i = -1
histos = []
snHistos = []
all1Histos = []
all8Histos = []
xed = []
yed = []
for galID in galID_list:
    i+=1
    for ion in ion_list:
        print galID, '\t', ion

        # Open the data 
        abs_file = galID+'.'+ion+'.bulk_abscells.dat'
        lognH, logT = np.loadtxt(file_loc+abs_file, skiprows=1, usecols=(7, 8), unpack=True)

        # Bin the data
        H, xedges, yedges = np.histogram2d( lognH, logT, bins=numbins)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        
        # Rotate and filp the histogram
        H = np.rot90(H)
        H = np.flipud(H)        
        xed.append(xedges)
        yed.append(yedges)
        histos.append(H)
#        print H.shape() 
        if feedback_list[i]=='dwSN':
            snHistos.append(H)
        elif feedback_list[i]=='dwALL_1':
            all1Histos.append(H)
        elif feedback_list[i]=='dwALL_8':
            all8Histos.append(H)
        else:
            print 'Unknown feedback name: ', feedback_list[i]
            sys.exit()


# Normalize dwALL_1 and dwALL_8 to the dwSN counts if
# the normalize flag is high
if normalize==1:
    for i in range(0,len(all1Histos)):
        allH = all1Histos[i]
        snH = snHistos[i]
        print allH
        print allH.shape
        print ''
        print snH
        print snH.shape
        print ''
        print (snH==0).shape
        with np.errstate(divide='ignore'):
            allH = allH/snH
            allH[ snH == 0 ] = 0.
        all1Histos[i] = allH
    for i in range(0,len(all8Histos)):
        allH = all8Histos[i]
        snH = snHistos[i]
        with np.errstate(divide='ignore'):
            allH = allH/snH
            allH[snH==0] = 0.
        all8Histos[i] = allH
    
# Create the subplots
fig,((p11,p12,p13), (p21,p22,p23), (p31,p32,p33), (p41,p42,p43)) = plt.subplots(4,3,figsize=(10.2,10.8))
plot_list = [p11, p21, p31, p41, p12, p22, p32, p42, p13, p23, p33, p43]

labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']

labels = ['(a)', '(d)', '(g)', '(j)', 
          '(b)', '(e)', '(h)', '(k)', 
          '(c)', '(f)', '(i)', '(l)']
for i in range(0,len(plot_list)):
    ax = plot_list[i]
    H = histos[i]
#            Hmasked = np.ma.masked_where(H==0,H)
    if i<4:     
        # Sn galaxy
        H = snHistos[i]
    elif i>=4 and i<8:
        # All1 galaxy
        H = all1Histos[i-4]
    else:
        # ALL8 galaxy
        H = all8Histos[i-8]

    if i==0 or i==4 or i==8:
        # Ion is HI
        cbarmin = 0
        cbarmax = 4
        cbarLabel = '$\log$ (Counts)'
    elif i==1 or i==5 or i==9:
        # Ion is MgII
        cbarmin = 0
        cbarmax = 3
        cbarLabel = '$\log$ (Counts)'
    elif i==2 or i==6 or i==10:
        # Ion is CIV
        cbarmin = 0
        cbarmax = 3.5
        cbarLabel = '$\log$ (Counts)'
    elif i==3 or i==7 or i==11:
        cbarmin = 0
        cbarmax = 3.5
        cbarLabel = '$\log$ (Counts))'

    H = np.ma.masked_where(H==0,H)
    H = np.log10(H)
    xedges = xed[i]
    yedges = yed[i]
#    mesh = ax.pcolormesh(xedges, yedges, H, vmin=cbarmin, vmax=cbarmax)
    mesh = ax.pcolormesh(xedges, yedges, H)
    ax.set_xlim([-8, 1])
    ax.set_ylim([2,8])
    ax.text(-1, 7, labels[i])
#    xlabels = ax.get_xticklabels()
#    ax.set_xticklabels(xlabels, fontsize=14)
    if i==0:
        ax.set_title('dwSN', size=12)
        ax.set_ylabel('HI \n $\log$ (T) [K] ', multialignment='center')
    elif i==1:
        ax.set_ylabel('MgII \n $\log$ (T) [K] ', multialignment='center')
    elif i==2:
        ax.set_ylabel('CIV \n $\log$ (T) [K] ', multialignment='center')
    elif i==3:
        ax.set_ylabel('OVI \n $\log$ (T) [K] ', multialignment='center')
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==4:
        ax.set_title('dwALL_1', size=12)
    elif i==7:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==8:
        ax.set_title('dwALL_8', size=12)
    elif i==11:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')

    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(cbarLabel, rotation=270, fontsize=12)
    
dropx = [p12, p22, p32, p13, p23, p33, p11, p21, p31]
plt.setp([a.get_xticklabels() for a in dropx],visible=False)

dropy = [p12, p22, p32, p13, p23, p33, p42, p43]
plt.setp([a.get_yticklabels() for a in dropy],visible=False)

plt.tight_layout()
if normalize==0:
    s = file_loc+'abscell_phase.master_bulk{0:d}.pdf'.format(numbins)
else:
    s = file_loc+'abscell_phase.master_bulk{0:d}_normed.pdf'.format(numbins)

plt.savefig(s, bbox_inches='tight')
