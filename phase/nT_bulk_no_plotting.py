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
from compHist import *

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

print '\nStart Binning...'


i = -1
histos = []
xed = []
yed = []
for galID in galID_list:
    i+=1
    
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
#        Hmasked = np.ma.masked_where(H==0,H)
        
        # Take the log of the count
#        H = np.log10(H)
        
        xed.append(xedges)
        yed.append(yedges)
        histos.append(H)
        
        # Write the histogram to file
        filename = file_loc+galID+'.'+ion+'.abscell.hist'
        f = open(filename, 'w')
        for k in range(0,len(H[:,0])):
            s = ''
            for l in range(0,len(H[0:,0])):
                s += '{0:.5f}\n'.format(H[k,l])
            f.write(s)

        f.close()


# Determine the chi2 goodness-of-fit for each ion
# Take the SN-only run (o2) as the 'expected'
# Ion    Indicies in histos
# HI     0, 4, 8 
# MgII   1, 5, 9
# CIV    2, 6, 10
# OVI    3, 5, 11

print '\n'
for i in range(0,len(ion_list)):

    snHist = histos[i]
    a1Hist = histos[i+4]
    a8Hist = histos[i+8]

    a1Chi2 = chi2GOF( snHist, a1Hist)
    a8Chi2 = chi2GOF( snHist, a8Hist)

    a1Chi2mat = chi2Matrix( snHist, a1Hist)
    a8Chi2mat = chi2Matrix( snHist, a8Hist)

    a1Chstwo, a1prob, a1df = chstwo( snHist, a1Hist)
    a8Chstwo, a8prob, a8df = chstwo( snHist, a8Hist)


    print '\nIon: ', ion_list[i]
#    print 'ALL_1 Chi2 GOF   : ', a1Chi2
#    print 'ALL_8 Chi2 GOF   : ', a8Chi2
#    print 'ALL_1 Chi2 Matrix: ', a1Chi2mat
#    print 'ALL_8 Chi2 Matrix: ', a8Chi2mat
    print 'ALL_1 Chstwo: ', a1Chstwo, a1prob, a1df
    print 'ALL_8 Chstwo: ', a8Chstwo, a8prob, a8df


    
