#!/usr/bin/python
#
# Filename: quick_cover.py
# Version: 1
#
# Author: Jacob Vander Vliet
# Date: 20/02/2014
#
# Compute the non-cumulative covering fraction
# Cut the impact space into annulli
# Add a loop over all EW cuts
# Add in errors
#
# Plots the covering fraction for a
# single ion for a single galaxy
#
# Usage:
#  python coverfrac.py <ALL file> <Rvir>  <ew cut>


# Function for determining vertical error bars
# Uses incomplete Beta functions
def funcup(x,n1,n2,cl):
    a = n1 + 1.0
    b = n2
    return cl - sc.betainc(a,b,x)

def funcdn(x,n1,n2,cl):
    a = n2 + 1.0
    b = n1
    return cl - sc.betainc(a,b,x)


import numpy as np
import scipy.special as sc
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import math
import matplotlib
import os
import subprocess as sp


font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font',**font)

allfile = sys.argv[1]
rvir = float(sys.argv[2])
ewcut = float(sys.argv[3])

galID = allfile.split('.')[0]
ion = allfile.split('.')[1]

# Define the radial bins to use
nbins = 6
binsize = 1.5/nbins


# Determine the bin edges
Dmax = []
Dmin = [0.0]
for j in range(0,nbins):
    Dmax.append( (j+1)*binsize )
    Dmin.append( (j+1)*binsize )
del Dmin[-1]

los=[]
CL = 0.8413
tol = 1e-4
pmin = 0.0
pmax = 1.0
iter = 1000

# Open file with list of all LOS
f = open('lines.info','r')
f.readline()
f.readline()
for line in f:
    los.append(line.split())
f.close()

f = open(allfile,'r')
header = (f.readline()).split()
    
impact=[]
abs=[]
absimpact=[]
covering=[]
imp=[]
horizerrneg=[]
horizerrpos=[]
verterrneg=[]
verterrpos=[]
for line in f:
    abs.append(line.split())
f.close()

# Loop over the different impact parameters
for j in range(0,len(Dmax)):
    maxrad = Dmax[j]
    minrad = Dmin[j]
    hit = 0.0
    total = 0.0
    xaxispoint = 0.0

    # Loop over all lines from the .sysabs file to get the
    # number of lines with significant absorption (SL > 3)
    # with an impact parameter between minrad and maxrad
    for k in range(0,len(abs)):
        absimpact.append(float(abs[k][1]) / rvir)
        width = float(abs[k][5])
        if absimpact[k]>minrad and absimpact[k]<maxrad and width>ewcut:
            hit += 1
            xaxispoint += absimpact[k]
            
    # Loop over all lines from impactparam.dat to get the
    # total number of sightlines with an impact parameter
    # between minrad and maxrad
    for k in range(0,len(los)):
        impact.append(float(los[k][1]) / rvir)
        if impact[k] > minrad and impact[k] < maxrad:
            total += 1
                    
    fraction = hit/total
            
    covering.append(fraction)
                    
    # Determine the x location of the point. If hit>0, then the location is the
    # average impact parameter in that bin. If hit=0, then the location is the
    # midpoint of the bin.
    if hit > 0.0:
        imp.append(xaxispoint/hit)  # location of point on x-axis
    else:
        imp.append((maxrad-minrad)/2.0 + minrad) 
            
    # Determine error bars
    #   Horizontal error bars are width of bin
    horizerrneg.append(imp[j]-minrad)
    horizerrpos.append(maxrad-imp[j])

    #   Vertical error bars are found using the incomplete beta function
    top = opt.brentq(funcup,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)
    bot = 1.0 - opt.brentq(funcdn,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)

    if fraction == 1.0:
        top = 1.0
    if fraction == 0.0:
        bot = 0.0

    verterrpos.append(top-fraction)
    verterrneg.append(fraction-bot)


plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos],marker='v',c='k',linestyle='none')
    
#plt.title(galID +'; '+ion+ '; EW Cut: {0:0.1f} $\AA$'.format(ewcut) )
plt.title(ion+ '; EW Cut: {0:0.1f} $\AA$'.format(ewcut) )
plt.xlabel('Impact Parameter (Rvir)')
plt.ylabel('Covering Fraction')
plt.ylim([-0.1,1.1])
plt.xlim([min(Dmin),max(Dmax)])
plt.xscale('log')
#plt.legend(numpoints=1,loc=1,prop={'size':10}, frameon=False)
plt.gca().get_frame().set_linewidth(2) 



plt.tight_layout()
plt.draw()
plt.savefig(galID+'_'+ion+'_{0:0.1f}_cover_quick.pdf'.format(ewcut))
