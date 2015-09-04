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
gasfile_loc = '/home/jacob/matrix/sebass_gals/dwarfs/gasfiles/'

# Binning settings
numbins = 50
mindense = -8.0
maxdense = 1.0
mintemp = 2.0
maxtemp = 8.0

# Create bins
n_bins = np.linspace(mindense, maxdense, num=numbins)
t_bins = np.linspace(mintemp, maxtemp, num=numbins)


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

        # Create arrays to hold the ionization counts
        total_count = np.zeros((numbins, numbins))
        coll_count = np.zeros((numbins, numbins))
        phot_count = np.zeros((numbins, numbins))
        eqal_count = np.zeros((numbins, numbins))

        # Get the index of each cell in the bins
        n_index = np.digitize(lognH, n_bins, right=True)
        t_index = np.digitize(logT,  t_bins, right=True)

        print n_index.max()

        # Clip the index arrays to not include indicies beyond 
        # the scope of the array
        np.clip(n_index, 0, numbins-1, out=n_index)
        np.clip(t_index, 0, numbins-1, out=t_index)

        print n_index.max()


        # Loop over the cells to determine the fraction of cells in
        # each bin that is collisionally ionized
        for j in range(0,len(lognH)):

            dense = lognH[j]
            temp = logT[j]
            photo_time = tph[j]
            coll_time = tcoll[j]
            n_ind = n_index[j]
            t_ind = t_index[j]

            if round(dense)==-4 and round(temp)==5:
                print dense, temp, n_ind, t_ind, photo_time, coll_time, coll_time<photo_time

            # Determine what is the dominant ionization mechanism
            if coll_time<photo_time:
                coll_count[n_ind, t_ind] += 1.0
            elif coll_time>photo_time:
                phot_count[n_ind, t_ind] += 1.0
            else:
                eqal_count[n_ind, t_ind] += 1.0
            total_count[n_ind, t_ind] += 1.0

        
            
        print '\nColl_count:'
        print coll_count.max()
        print '\nTotal_count:'
        print total_count.max()

        # Determine the fraction of cells in each bin that 
        # are collisionally ionized
        coll_fraction = np.divide(coll_count, total_count)
        coll_fraction = np.nan_to_num(coll_fraction)

        phot_fraction = np.divide(phot_count, total_count)
        phot_fraction = np.nan_to_num(phot_fraction)

        eqal_fraction = np.divide(eqal_count, total_count)
        eqal_fraction = np.nan_to_num(eqal_fraction)

        print '\nColl_Fraction:'
        print coll_fraction.min()
        print coll_fraction.max()
        print coll_fraction.mean()



        print '\nColl_Fraction:'
        print coll_fraction

        coll_fraction = np.ma.masked_where(coll_fraction==0, coll_fraction)
        coll_fraction = np.rot90(coll_fraction)
        coll_fraction = np.flipud(coll_fraction)
        
        phot_fraction = np.ma.masked_where(phot_fraction==0, phot_fraction)
        phot_fraction = np.rot90(phot_fraction)
        phot_fraction = np.flipud(phot_fraction)

        eqal_fraction = np.ma.masked_where(eqal_fraction==0, eqal_fraction)
        eqal_fraction = np.rot90(eqal_fraction)
        eqal_fraction = np.flipud(eqal_fraction)
        
        
        # Plot (this is only for testing)
        fig = plt.figure()
        plt.pcolormesh(n_bins, t_bins, coll_fraction)
        plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
        plt.ylabel(' $\log$ (T) [K] ')
        plt.xlim([-8, 1])
        plt.ylim([2,8])
        plt.title(feedback_list[i] +', '+ ion)
        cbar = plt.colorbar()
#        cbar.ax.set_ylabel('$\log$ (Counts)')

        plt.savefig('dum'+ion+'coll.pdf')

        plt.cla()
        plt.clf()

        fig = plt.figure()
        plt.pcolormesh(n_bins, t_bins, phot_fraction)
        plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
        plt.ylabel(' $\log$ (T) [K] ')
        plt.xlim([-8, 1])
        plt.ylim([2,8])
        plt.title(feedback_list[i] +', '+ ion)
        cbar = plt.colorbar()
#        cbar.ax.set_ylabel('$\log$ (Counts)')

        plt.savefig('dum'+ion+'phot.pdf')

        plt.cla()
        plt.clf()

     
        fig = plt.figure()
        plt.pcolormesh(n_bins, t_bins, eqal_fraction)
        plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
        plt.ylabel(' $\log$ (T) [K] ')
        plt.xlim([-8, 1])
        plt.ylim([2,8])
        plt.title(feedback_list[i] +', '+ ion)
        cbar = plt.colorbar()
#        cbar.ax.set_ylabel('$\log$ (Counts)')

        plt.savefig('dum'+ion+'eqal.pdf')

        plt.cla()
        plt.clf()

     
    sys.exit()
