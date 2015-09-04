#!/usr/bin/python
# Filename: nT.py
#
# Author: Jacob Vander Vliet
# Version: 2
# Date: 4/12/2013
# Last Modified: 4/12/13
#
# Description:
#   Creates the phase diagram of the gas in the general box (before rates)
#
# Instructions:
#   python nT.py <name of small box> <inner radius> <outer radius>
#
# For everything inside 5 kpc:
#   python nT.py dwarf9o_a1.001.txt 0 5
#
# For everythin outside of 5 kpc:
#   python nT.py dwarf9o_a1.001.txt 5 inf
#
# If, instead, you want to cut the data a percent of Rvir
#   python nT.py <percent>


import numpy as np
import matplotlib.pyplot as plt
from os import system
import scipy.interpolate
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib as mpl
import sys
import subprocess as sp
import os


verbose = 0

print len(sys.argv)


sp.call('ls *SZ* > files.list', shell=True)
f = open('files.list')

exp = []
for line in f:
    l = line.split('a')[2]
    str1 = l.split('.')[0]
    str2 = l.split('.')[1]
    exp.append(str1+'.'+str2)

f.close()

# Go into each directory
for i in range(0,len(exp)):


    os.chdir('./a'+exp[i])
    
    # Get the name of the galaxy box
    com = 'ls *GZ*.txt > name.list'
    sp.call(com, shell=True)
    f = open('name.list')
    filename = f.readline().strip()
    f.close()

    # Need to get Rvir
    sumloc = '/home/matrix3/jrvander/galaxy_files/summaries/'
    # Get galID and expansion parameter from filename
    galID = filename.split('_')[0]
    sumloc = sumloc+galID+'.dat'
#        str1 = sys.argv[1].split('_')[1]
#        str2 = (str1.split('.')[0]).split('a')[1]
#        str3 = str1.split('.')[1]
#        aexpn = str2+'.'+str3
        
    f = open(sumloc)
    f.readline()
    f.readline()

#        print aexpn
    found=0
    for line in f:
        if found==0:
            a = line.split()
                
            if float(a[0])==float(exp[i]):
                found=1
                rvir = float(line.split()[3])
    f.close()
    cutoff = float(sys.argv[1])/100.*rvir


#    filename = sys.argv[1]
    galname = filename.split('_')[0]

    data = np.loadtxt(filename,skiprows=2)

    # Pull out the cell location
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]

    # Pull out density and temp and take their log
    density = data[:,7]
    temp = data[:,8]
    
    logn = [math.log10(k) for k in density]
    logT = [math.log10(k) for k in temp]


    # Select the cells that are within the right range
    nlow = []
    Tlow = []
    nhigh = []
    Thigh = []
    for j in range(0,len(x)):
        dist = math.sqrt( x[j]*x[j] + y[j]*y[j] + z[j]*z[j] )

#    if dist<outer and dist>inner:
#        n.append(logn[i])
#        T.append(logT[i])
        if dist<cutoff:
            nlow.append(logn[j])
            Tlow.append(logT[j])
        else:
            nhigh.append(logn[j])
            Thigh.append(logT[j])
    


    if verbose==1:
        print 'Min density: ',density.min()
        print 'Max density: ',density.max()
        print 'Min temp:    ',temp.min()
        print 'Max temp:    ',temp.max() 

        print ''
        print 'Min logn: ',min(logn)
        print 'Max logn: ',max(logn)
        print 'Min logT: ',min(logT)
        print 'Max logT: ',max(logT)
        
        print ''
        print 'Temp at min logn: ',logT[logn.index(min(logn))]
        print 'Temp at max logn: ',logT[logn.index(max(logn))]


    # Bin the data
    numbins = 50
    Hl, xedgesl, yedgesl = np.histogram2d( nlow, Tlow, bins=numbins)
    Hh, xedgesh, yedgesh = np.histogram2d( nhigh, Thigh, bins=numbins)
    extentl = [xedgesl[0], xedgesl[-1], yedgesl[0], yedgesl[-1]]
    extenth = [xedgesh[0], xedgesh[-1], yedgesh[0], yedgesh[-1]]
    
    # Rotate and filp the histogram
    Hl = np.rot90(Hl)
    Hl = np.flipud(Hl)
    Hh = np.rot90(Hh)
    Hh = np.flipud(Hh)

    # Mask the bins where the count is zero
    Hmaskedl = np.ma.masked_where(Hl==0,Hl)
    Hmaskedh = np.ma.masked_where(Hh==0,Hh)
    
    # Take the log of the count
    Hmaskedl = np.log10(Hmaskedl)
    Hmaskedh = np.log10(Hmaskedh)


    # Plot and save the figure
    fig = plt.figure()
    plt.pcolormesh(xedgesl,yedgesl,Hmaskedl)
    plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    plt.ylabel(' $\log$ (T) [K] ')
    plt.xlim([-8, 1])
    plt.ylim([2,8])
    plt.title('Inside {0}% Rvir ({1:0.2f} kpc)'.format(sys.argv[1],cutoff))
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('$\log$ (Counts)')
    fig.savefig('phase_'+galname +'_'+exp[i]+'_in_{0}percent.pdf'.format(sys.argv[1]))
    
    
    plt.cla()
    plt.clf()
    
    fig = plt.figure()
    plt.pcolormesh(xedgesh,yedgesh,Hmaskedh)
    plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    plt.ylabel(' $\log$ (T) [K] ')
    plt.xlim([-8, 1])
    plt.ylim([2,8])
    plt.title('Outside {0}% Rvir ({1:0.2f} kpc)'.format(sys.argv[1],cutoff))
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('$\log$ (Counts)')
    fig.savefig('phase_'+galname +'_'+exp[i]+'_out_{0}percent.pdf'.format(sys.argv[1]))

    
    os.chdir('..')
