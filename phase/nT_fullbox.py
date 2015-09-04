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

gasfile = sys.argv[1]
galID = gasfile.split('_')[0]
e1 = gasfile.split('.')[0].split('a')[1]
e2 = gasfile.split('.')[1]
expn = e1+'.'+e2

data = np.loadtxt(gasfile,skiprows=2)

# Pull out the cell location
x = data[:,1]
y = data[:,2]
z = data[:,3]

# Pull out density and temp and take their log
density = data[:,7]
temp = data[:,8]

logn = [math.log10(k) for k in density]
logT = [math.log10(k) for k in temp]

print 'Temperature range: ', min(logT), max(logT)


# Bin the data
numbins = 50
H, xedges, yedges = np.histogram2d( logn, logT, bins=numbins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

# Rotate and filp the histogram
H= np.rot90(H)
H = np.flipud(H)

# Mask the bins where the count is zero
Hmasked = np.ma.masked_where(H==0,H)

# Take the log of the count
Hmasked = np.log10(Hmasked)

# Plot and save the figure
fig = plt.figure()
plt.pcolormesh(xedges,yedges,Hmasked)
plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
plt.ylabel(' $\log$ (T) [K] ')
plt.xlim([-8, 1])


plt.ylim([2,8])
cbar = plt.colorbar()
cbar.ax.set_ylabel('$\log$ (Counts)')
fig.savefig(galID+'.'+expn+'.fullboxphase.pdf')    

plt.cla()
plt.clf()

