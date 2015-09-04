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
#   python nT_bulk.py <galID> <expn> <Rvir>
#

import numpy as np
import matplotlib.pyplot as plt
import math
import sys



if len(sys.argv)!=4:
    print 'Usage:\n\t python nT_bulk_plotting.py <galID> <expn> <Rvir> \n'
    sys.exit()
galID = sys.argv[1]
expn = sys.argv[2]
rvir = float(sys.argv[3])

verbose = 0
bulk = 1

ion_list = ['HI', 'MgII', 'CIV', 'OVI']


numbins = 50
i = -1
histos = []
xed = []
yed = []
i+=1
for ion in ion_list:
    print galID, '\t', ion

    # Open the data 
    abs_file = './'+ion+'/'+galID+'.'+expn+'.'+ion+'.abs_cells.dat'
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

    
# Create the subplots
fig,((p11,p12), (p21,p22)) = plt.subplots(2,2,figsize=(10.2,10.2))
plot_list = [p11, p21, p12, p22]
labels = ['(a)', '(b)', '(c)', '(d)'] 
labels = ['(a)', '(c)', '(b)', '(d)'] 

for i in range(0,len(plot_list)):
    ax = plot_list[i]
    H = histos[i]
    if i==0:
        # Ion is HI
        cbarmin = 0
        cbarmax = 4
        cbarLabel = '$\log$ (Counts)'
    elif i==2:
        # Ion is MgII
        cbarmin = 0
        cbarmax = 3
        cbarLabel = '$\log$ (Counts)'
    else:
        # Ion is CIV or OVI
        cbarmin = 0
        cbarmax = 3.5
        cbarLabel = '$\log$ (Counts)'

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
        ax.set_ylabel('HI \n $\log$ (T) [K] ', multialignment='center')
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==2:
        ax.set_ylabel('MgII \n $\log$ (T) [K] ', multialignment='center')
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==1:
        ax.set_ylabel('CIV \n $\log$ (T) [K] ', multialignment='center')
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==3:
        ax.set_ylabel('OVI \n $\log$ (T) [K] ', multialignment='center')
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    
    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(cbarLabel, rotation=270, fontsize=12)
    
#dropx = [p12, p22, p32, p13, p23, p33, p11, p21, p31]
#plt.setp([a.get_xticklabels() for a in dropx],visible=False)

#dropy = [p12, p22, p32, p13, p23, p33, p42, p43]
#plt.setp([a.get_yticklabels() for a in dropy],visible=False)

plt.tight_layout()
s = file_loc+'abscell_phase_{0:s}_{1:s}_{2:d}.pdf'.format(galID, expn, numbins)
plt.savefig(s, bbox_inches='tight')
