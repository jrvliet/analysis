#!/usr/bin/python
#
# Filename: coverfrac_random.py
# Version: 1
#
# Author: Jacob Vander Vliet
# Date: 27/07/2014
#
# Compute the non-cumulative covering fraction
# Cut the impact space into annulli
# Run on the large ALL files stored in:
#  /home/matrix3/jrvander/sebass_gals/dwarfs/ALLfiles/
#
# Selects n random sightlines from the lines.info files
# Calcualtes the covering fraction for all ions for the
# Chosen lines
#
# Usage:
#  python coverfrac_random.py <ewcut> <number of lines>



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
import os.path as op
import random as rn


ion_list = ['HI', 'MgII', 'CIV', 'OVI']
galID_list = ['D9o2', 'D9q', 'D9m4a']
mark=['v','o','x']
colors=['fuchsia', 'blue', 'black']
line=['-','.']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]

los_per_box = 999
CL = 0.8413
tol = 1e-4
pmin = 0.0
pmax = 1.0
iter = 1000

# Define the EW cuts to use
ewcut = float(sys.argv[1])

# Set the number of lines of sight to use
nLOS = int(sys.argv[2])

# Define the radial bins to use
nbins = 15
binsize = 1.5/nbins

# Determine the bin edges
Dmax = []
Dmin = [0.0]
for j in range(0,nbins):
    Dmax.append( (j+1)*binsize )
    Dmin.append( (j+1)*binsize )
del Dmin[-1]


# Determine the random LOS to use by 
# generating n random numbers from 1-6000
rn.seed(25655)
total_losnum = los_per_box*len(expn_list_1)
skewers = []  # array of LOS nums
for i in range(0,n):
    skewers.append(randint(1,total_losnum))


for ion in ion_list:

    # Loop over the different feedback models
    for galID in galID_list:

        for los in skewers:

            # Determine which snapshot this is
            snapshot = los % los_per_box
            expansion = expn_list[galID_list.index(galID), snapshot]
        
            # Read in the lines file
            lines = np.loadtxt('lines.'+galID+'.a'+expansion+'.info', skiprows=2)
        
            # Get the info for this LOS
            los_imp = lines[los,1]
        
            # Loop over the ions

            
            # Read in the appropriate ALL file
            allf = np.loadtxt(galID+'.a'+expansion+'.'+ion+'.ALL.sysabs.large',skiprows=1)

            # Loop through the ALL file and determine if the line is in the ALL file
            found=0
            for j in range(0,len(allf[:,0])):
                if allf[j,0]==los:
                    found=1
            if found==1:
                hit+=1
                total+=1
            else:
                total+=1

            

for ion in ion_list:
    
    for galID in galID_list:
        allfile = galID+'.'+ion+'.ALL.sysabs.master'

        impact=[]
        abs=[]
        absimpact=[]
        covering=[]
        imp=[]
        horizerrneg=[]
        horizerrpos=[]
        verterrneg=[]
        verterrpos=[]

        # Read in the allfile
        f = open(allfile, 'r')
        f.readline()
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
            zerobinabs = []
            zerobin = []

            # Loop over all lines from the .sysabs file to get the
            # number of lines with significant absorption (SL > 3)
            # with an impact parameter between minrad and maxrad
            for k in range(0,len(abs)):
#                absimpact.append(float(abs[k][1]) / rvir)
                absimpact.append(float(abs[k][22]))
                width = float(abs[k][5])
                if absimpact[k]>minrad and absimpact[k]<maxrad and width>ewcut:
                    hit += 1
                    xaxispoint += absimpact[k]
                    if minrad==0.0:
                        zerobinabs.append( abs[k][0] )
            
            
            # Loop over all lines from impactparam.dat to get the
            # total number of sightlines with an impact parameter
            # between minrad and maxrad
            for k in range(0,len(los)):
                impact.append(float(los[k][1]) / rvir)
                if impact[k] > minrad and impact[k] < maxrad:
                    total += 1
                    if minrad==0.0:
                        continue
                    
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
            elif fraction == 0.0:
                bot = 0.0

            verterrpos.append(top-fraction)
            verterrneg.append(fraction-bot)

        plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos],label=galID,marker=mark[i],c=colors[i],linestyle='none')
        

    plt.xlabel('$\\rho$ / R$_{vir}$')
    ylab = 'C$_f$ ['+ion+']'
    plt.ylabel(ylab)
    plt.ylim([-0.1,1.1])
    plt.xlim([min(Dmin),1.0])
    plt.legend(numpoints=1,loc=1,prop={'size':15}, frameon=False)
    plt.gca().get_frame().set_linewidth(2) 
    
    plt.savefig('master_'+ion+'_{0:0.1f}mA_cover.eps'.format(ewcut*1000))



