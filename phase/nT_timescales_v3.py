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
feedback_list = ['dwSN', 'dwALL\_1', 'dwALL\_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']
#ion_list = ['OVI']


gasfile_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/gasfiles/'
gasfile_loc = '/home/jacob/matrix/sebass_gals/dwarfs/gasfiles/'
gasfile_loc = '/home/jacob/research/dwarf_gasfiles/'

abscell_loc = gasfile_loc.replace('gasfiles','abscells')


# Binning settings
numbins = 50
mindense = -8.0
maxdense = 1.0
mintemp = 2.0
maxtemp = 8.0

# Create bins
n_bins = np.linspace(mindense, maxdense, num=numbins)
t_bins = np.linspace(mintemp, maxtemp, num=numbins)

nstep = n_bins[0]-n_bins[1]
tstep = t_bins[0]-t_bins[1]

for i in range(0,len(galID_list)):

    galID = galID_list[i]

    print ''
    print galID
    # Create the four subplots
    fig,((p11,p12,p13), (p21,p22,p23), (p31,p32,p33), (p41,p42,p43)) = plt.subplots(4,3,figsize=(10.2,10.8))
    plot_list = [p11, p21, p31, p41, p12, p22, p32, p42, p13, p23, p33, p43]

    
    # Loop of over ions
    for ion in ion_list:
        print '\t', ion

        # Create arrays to hold the ionization counts
        total_count = np.zeros((numbins, numbins))
        coll_count = np.zeros((numbins, numbins))
        phot_count = np.zeros((numbins, numbins))
        eqal_count = np.zeros((numbins, numbins))
        alt_count = np.zeros((numbins, numbins))

        
        # Open the abs cells
        absfile = galID+'.'+ion+'.bulk_abscells.dat'
        absf = open(abscell_loc+absfile)
        absf.readline()
        
        old_expn = '0.0'
        for line in absf:

            expn = line.split()[13]

            # If this is different than the previous line, then 
            # a new  box needs to be read in
            if expn!=old_expn:
                # Read new box
                print '\t\tMoving to a = ',expn
                gasfile = galID+'_GZa'+expn+'.'+ion+'.txt'
                nH, T, tph, trec, tcoll = np.loadtxt(gasfile_loc+gasfile, skiprows=2, usecols=(7, 8, 17, 18, 19), unpack=True)
                
            old_expn = expn

            cellID = int(line.split()[2])
#            print cellID, len(nH)

            # Get the properties of this cell
            dense = np.log10(nH[cellID-1])
            temp = np.log10(T[cellID-1])
            photo_time = tph[cellID-1]
            coll_time = tcoll[cellID-1]
            rec_time = trec[cellID-1]

            # Calculate the alternative collisional timescale
            alt_coll = (rec_time / coll_time) / (1 + (rec_time / coll_time))

            # Get the index of this cell in the n, t bins
            n_index = floor((dense - mindense) / nstep)
            t_index = floor((temp - mintemp) / tstep)

            # Clip the index arrays to not include indicies beyond 
            # the scope of the array
            if n_index==numbins:
                n_index = numbins-1
            if t_index==numbins:
                t_index = numbins-1

            # Determine the dominate ionization mechanism for this cell
            if coll_time<photo_time:
                coll_count[n_index, t_index] += 1.0
            elif coll_time>photo_time:
                phot_count[n_index, t_index] += 1.0
            elif alt_coll<photo_time:
                alt_count[n_index, t_index] += 1.0
            else:
                eqal_count[n_index, t_index] += 1.0
            total_count[n_index, t_index] += 1.0

        
            
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

        alt_fraction = np.divide(alt_count, total_count)
        alt_fraction = np.nan_to_num(alt_fraction)

        for index, value in np.ndenumerate(total_count):

            if value>0 and coll_fraction[index]==0:
                coll_fraction[index] = 1.0e-5

            if value>0 and alt_fraction[index]==0:
                alt_fraction[index] = 1.0e-5


        coll_fraction = np.ma.masked_where(coll_fraction==0, coll_fraction)
        coll_fraction = np.rot90(coll_fraction)
        coll_fraction = np.fliplr(coll_fraction)
        coll_fraction = np.log10(coll_fraction)

        alt_fraction = np.ma.masked_where(alt_fraction==0, alt_fraction)
        alt_fraction = np.rot90(alt_fraction)
        alt_fraction = np.fliplr(alt_fraction)
        alt_fraction = np.log10(alt_fraction)
        
#        phot_fraction = np.ma.masked_where(phot_fraction==0, phot_fraction)
#        phot_fraction = np.rot90(phot_fraction)
#        phot_fraction = np.flipud(phot_fraction)

#        eqal_fraction = np.ma.masked_where(eqal_fraction==0, eqal_fraction)
#        eqal_fraction = np.rot90(eqal_fraction)
#        eqal_fraction = np.flipud(eqal_fraction)
        
        
        # Plot (this is only for testing)
        fig = plt.figure()
        plt.pcolormesh(n_bins, t_bins, coll_fraction, vmin=-2, vmax=0)
        plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
        plt.ylabel(' $\log$ (T) [K] ')
        plt.xlim([-8, 1])
        plt.ylim([2,8])
        plt.title(feedback_list[i] +', '+ ion)
        cbar = plt.colorbar()
#        cbar.ax.set_ylabel('$\log$ (Counts)')
        plt.savefig(galID+'_'+ion+'_coll_{0:d}.pdf'.format(numbins))
        plt.cla()
        plt.clf()

        fig = plt.figure()
        plt.pcolormesh(n_bins, t_bins, alt_fraction, vmin=-2, vmax=0)
        plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
        plt.ylabel(' $\log$ (T) [K] ')
        plt.xlim([-8, 1])
        plt.ylim([2,8])
        plt.title(feedback_list[i] +', '+ ion)
        cbar = plt.colorbar()
#        cbar.ax.set_ylabel('$\log$ (Counts)')
        plt.savefig(galID+'_'+ion+'_alt_{0:d}.pdf'.format(numbins))
        plt.cla()
        plt.clf()

     

