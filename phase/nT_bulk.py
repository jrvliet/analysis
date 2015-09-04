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


verbose = 0


normalize = int(sys.argv[1])

if len(sys.argv)!=3:
    bulk = 0
else:
    bulk = int(sys.argv[2])


galID_list = ['D9o2', 'D9q', 'D9m4a']
feedback_list = ['dwSN', 'dwALL_1', 'dwALL_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']


file_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/abscells/'

i = -1
histos = []
xed = []
yed = []
for galID in galID_list:
    i+=1
    
    if bulk!=0:
        # Create the four subplots
        fig,((p11,p12,p13), (p21,p22,p23), (p31,p32,p33), (p41,p42,p43)) = plt.subplots(4,3,figsize=(10.2,10.8))
#        plot_list = [p11,p12,p13, p21,p22,p23, p31,p32,p33, p41,p42,p43]
        plot_list = [p11, p21, p31, p41, p12, p22, p32, p42, p13, p23, p33, p43]

    for ion in ion_list:
        print galID, '\t', ion

        # Open the data 
        abs_file = galID+'.'+ion+'.bulk_abscells.dat'
        lognH, logT = np.loadtxt(file_loc+abs_file, skiprows=1, usecols=(7, 8), unpack=True)
        print lognH[0], logT[0]

        # Bin the data
        numbins = 50
        H, xedges, yedges = np.histogram2d( lognH, logT, bins=numbins, normed=normalize)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        
        # Rotate and filp the histogram
        H = np.rot90(H)
        H = np.flipud(H)
        
        # Mask the bins where the count is zero
        Hmasked = np.ma.masked_where(H==0,H)
        
        # Take the log of the count
        Hmasked = np.log10(Hmasked)
        
        # Plot and save the figure
        if bulk==0:
            fig = plt.figure()
            plt.pcolormesh(xedges,yedges,Hmasked)
            plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
            plt.ylabel(' $\log$ (T) [K] ')
            plt.xlim([-8, 1])
            plt.ylim([2,8])
            plt.title(feedback_list[i] +', '+ ion)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('$\log$ (Counts)')
            if normalize==0:
                fig.savefig(file_loc+'abscell_phase.'+galID+'.'+ion+'.pdf')
                fig.savefig(file_loc+'abscell_phase.'+galID+'.'+ion+'.eps')
                fig.savefig(file_loc+'abscell_phase.'+galID+'.'+ion+'.jpg')
            else:
                fig.savefig(file_loc+'abscell_phase.'+galID+'.'+ion+'.normalized.pdf')
                fig.savefig(file_loc+'abscell_phase.'+galID+'.'+ion+'.normalized.eps')
                fig.savefig(file_loc+'abscell_phase.'+galID+'.'+ion+'.normalized.jpg')

            plt.cla()
            plt.clf()

        else:
            xed.append(xedges)
            yed.append(yedges)
            histos.append(Hmasked)
            continue

#            subplotnum = 221+ion_list.index(ion)
            ax = plot_list[ion_list.index(ion)]
#            fig = plt.subplot(subplotnum)
#            fig = plt.figure()
            mesh = ax.pcolormesh(xedges,yedges,Hmasked)
            ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
            ax.set_ylabel(' $\log$ (T) [K] ')
            ax.set_xlim([-8, 1])
            ax.set_ylim([2,8])
            ax.text(-1, 7, ion)
#            plt.title(feedback_list[i] +', '+ ion)
            cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
            cbar.ax.set_label('$\log$ (Counts)')
        


#    if bulk!=0:
#        plt.subplots_adjust(wspace=0.35, hspace=0.2)
#        plt.subplots_adjust(right=0.85, left=0.1)
#        plt.tight_layout()
#        s = file_loc+'abscell_phase.'+galID+'.master_bulk.eps'
#        plt.savefig(s, bbox_inches='tight')

#        for j in range(0,4):
#            plt.subplot(221+j)
#            plt.cla()
#            plt.clf()


# Plot if bulk is high
labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']

labels = ['(a)', '(d)', '(g)', '(j)', 
          '(b)', '(e)', '(h)', '(k)', 
          '(c)', '(f)', '(i)', '(l)']
if bulk!=0:
    for i in range(0,len(plot_list)):
        ax = plot_list[i]
        H = histos[i]
        xedges = xed[i]
        yedges = yed[i]
        mesh = ax.pcolormesh(xedges,yedges,H)
        ax.set_xlim([-8, 1])
        ax.set_ylim([2,8])
        ax.text(-1, 7, labels[i])
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
        cbar.ax.set_label('$\log$ (Counts)')
    
    dropx = [p12, p22, p32, p13, p23, p33, p11, p21, p31]
    plt.setp([a.get_xticklabels() for a in dropx],visible=False)

    dropy = [p12, p22, p32, p13, p23, p33, p42, p43]
    plt.setp([a.get_yticklabels() for a in dropy],visible=False)

    # Label the rows
#    fig.text(0.05, 0.85, 'HI', rotation='vertical')
#    fig.text(0.05, 0.63, 'MgII', rotation='vertical')
#    fig.text(0.05, 0.39, 'CIV', rotation='vertical')
#    fig.text(0.05, 0.15, 'OVI', rotation='vertical')

    plt.tight_layout()
    s = file_loc+'abscell_phase.master_bulk.eps'
    plt.savefig(s, bbox_inches='tight')
