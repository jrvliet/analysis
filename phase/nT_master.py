#!/usr/bin/python
# Filename: nT_master.py
#
# Author: Jacob Vander Vliet
# Version: 1
# Date: 29/07/2013
# Last Modified: 29/07/13
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
#   python nT.py dwarf9o_a1.001.txt p <percent>


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

print len(sys.argv)

# Check which method to use
if sys.argv[2]!='p':
    # Read in the limits
    if len(sys.argv)==2:
        inner = 0.
        outer = 1.0e6
    else:
        inner = float(sys.argv[2])
        outer = float(sys.argv[3])
        if outer=='inf':
            outer = 1.0e6
        else:
            outer = float(outer)

    cutoff

else:
    # Need to get Rvir
    sumloc = '/home/matrix3/jrvander/galaxy_files/summaries/'
    # Get galID and expansion parameter from filename
    galID = sys.argv[1].split('_')[0]
    sumloc = sumloc+galID+'.dat'
    str1 = sys.argv[1].split('_')[1]
    str2 = (str1.split('.')[0]).split('a')[1]
    str3 = str1.split('.')[1]
    aexpn = str2+'.'+str3
    
    f = open(sumloc)
    f.readline()
    f.readline()

    print aexpn
    found=0
    for line in f:
        if found==0:
            a = line.split()
            
            if float(a[0])==float(aexpn):
                found=1
                rvir = float(line.split()[3])
    f.close()
    cutoff = float(sys.argv[3])/100.*rvir


#print 'Inner radius = ',inner
#print 'Outer radius = ',outer
# Get name of galaxy box from command arguement and read in data
filename = sys.argv[1]
galname = filename.split('_')[0]

data = np.loadtxt(filename,skiprows=2)

# Pull out the cell location
x = data[:,1]
y = data[:,2]
z = data[:,3]

# Pull out density and temp and take their log
density = data[:,7]
temp = data[:,8]

logn = [math.log10(i) for i in density]
logT = [math.log10(i) for i in temp]


# Select the cells that are within the right range
nlow = []
Tlow = []
nhigh = []
Thigh = []
for i in range(0,len(x)):
    dist = math.sqrt( x[i]*x[i] + y[i]*y[i] + z[i]*z[i] )

#    if dist<outer and dist>inner:
#        n.append(logn[i])
#        T.append(logT[i])
    if dist<cutoff:
        nlow.append(logn[i])
        Tlow.append(logT[i])
    else:
        nhigh.append(logn[i])
        Thigh.append(logT[i])
    


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
plt.title('Inside {0}% Rvir ({1:0.2f} kpc)'.format(sys.argv[3],cutoff))
cbar = plt.colorbar()
cbar.ax.set_ylabel('$\log$ (Counts)')
fig.savefig('phase_'+galname +'inside_{0}_{1}.pdf'.format(sys.argv[3],cutoff))


plt.cla()
plt.clf()

fig = plt.figure()
plt.pcolormesh(xedgesh,yedgesh,Hmaskedh)
plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
plt.ylabel(' $\log$ (T) [K] ')
plt.xlim([-8, 1])
plt.ylim([2,8])
plt.title('Outside {0}% Rvir ({1:0.2f} kpc)'.format(sys.argv[3],cutoff))
cbar = plt.colorbar()
cbar.ax.set_ylabel('$\log$ (Counts)')
fig.savefig('phase_'+galname +'outside_{0}_{1}.pdf'.format(sys.argv[3],cutoff))
