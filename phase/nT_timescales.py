#!/usr/bin/python
# Filename: nT_timescales.py
#
# Author: Jacob Vander Vliet
# Version: 1
# Date: 7/10/2014
# Last Modified: 7/10/2014
#
# Description:
#   Plots all cells in the box on the nT space
#   where the color corresponds the fraction of cells
#   in that bin that are collisionally ionized
#
# Instructions:
#   python nT_timescales.py
#

import numpy as np
import matplotlib.pyplot as plt
import sys
from math import floor

np.seterr(divide='ignore')

galID_list = ['D9o2', 'D9q', 'D9m4a']
feedback_list = ['dwSN', 'dwALL_1', 'dwALL_8']
#ion_list = ['HI', 'MgII', 'CIV', 'OVI']
ion_list = ['OVI']


gasfile_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/gasfiles/'
numbins = 50

for i in range(0,len(galID_list)):

    galID = galID_list[i]

    # Create the four subplots
    fig,((p11,p12,p13), (p21,p22,p23), (p31,p32,p33), (p41,p42,p43)) = plt.subplots(4,3,figsize=(10.2,10.8))
    plot_list = [p11, p21, p31, p41, p12, p22, p32, p42, p13, p23, p33, p43]

    
    # Loop of over ions
    for ion in ion_list:
        
        # Open the data
        gasfile = galID+'_'+ion+'_GZ_combined.txt'
        
        nH, T, tph, tcoll = np.loadtxt(gasfile_loc+gasfile, skiprows=2, usecols=(7, 8, 17, 19), unpack=True)
        lognH = np.log10(nH)
        logT = np.log10(T)

        # Bin the data
        H, xedges, yedges = np.histogram2d( lognH, logT, bins=numbins)

        dense_step = xedges[1]-xedges[0]
        temp_step = yedges[1]-yedges[0]
        dense_min = xedges[0]
        dense_max = xedges[-1]
        temp_min = yedges[0]
        temp_max = yedges[-1]


#        print 'Density:'
#        print xedges
#        print dense_step
#        print dense_min
#        print dense_max
#        print dense_min + 50*dense_step

        
        print len(H[:,0])
        print len(H[0,:])


        # Create arrays to hold the ionization counts
        coll_count = np.zeros_like(H)
        total_count = np.zeros_like(H)
        coll_fraction = np.zeros_like(H)

        total_count = np.zeros((numbins, numbins))
        coll_count = np.zeros((numbins, numbins))

        print H
        print coll_count
        print total_count


        # Loop over the cells 
        for j in range(0,len(lognH)):

            dense = lognH[j]
            temp = logT[j]

            n_index = floor((dense - dense_min) / dense_step)
            t_index = floor((temp - temp_min) / temp_step)
            
            # If the temp or density is at the maximum level, the 
            # derived indecies will be too high by one.
            if temp==temp_max:
                t_index = t_index-1
            if dense==dense_max:
                n_index = n_index-1


            # Determine what is the dominant ionization mechanism
            photo_time = tph[j]
            coll_time = tcoll[j]

            if coll_time<photo_time:
                coll_count[n_index, t_index] += 1.0
            total_count[n_index, t_index] += 1.0

            
        print '\nColl_count:'
        print coll_count
        print '\nTotal_count:'
        print total_count

        # Determine the fraction of cells in each bin that 
        # are collisionally ionized
        coll_fraction = np.divide(coll_count, total_count)

        print '\nColl_Fraction:'
        print coll_fraction

        coll_fraction = np.nan_to_num(coll_fraction)

        print '\nColl_Fraction:'
        print coll_fraction

        print xedges
        print yedges

        coll_fraction = np.ma.masked_where(coll_fraction==0, coll_fraction)
#        coll_fraction = np.rot90(coll_fraction)
#        coll_fraction = np.flipud(coll_fraction)
        
        
        # Plot (this is only for testing)
        fig = plt.figure()
        plt.pcolormesh(xedges,yedges,coll_fraction)
        plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
        plt.ylabel(' $\log$ (T) [K] ')
        plt.xlim([-8, 1])
        plt.ylim([2,8])
        plt.title(feedback_list[i] +', '+ ion)
        cbar = plt.colorbar()
#        cbar.ax.set_ylabel('$\log$ (Counts)')

        plt.savefig('dum'+ion+'.pdf')

        plt.cla()
        plt.clf()

     
    sys.exit()
