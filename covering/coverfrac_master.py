#!/usr/bin/python
#
# Filename: coverfrac_master.py
# Version: 1
#
# Author: Jacob Vander Vliet
# Date: 27/07/2014
#
# Compute the non-cumulative covering fraction
# Cut the impact space into annulli
# Add a loop over all EW cuts
# Add in errors
# Fixed the name of galaxies
# Works on the master ALL files that combine results across 
# multiple snapshots
#
# Usage:
#  python coverfrac_master.py <ewcut> <bulk?>
#
# Bulk flag: if=1, will plot all four ions on one plot



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

if len(sys.argv)!=3:
    bulk = 0
else:
    bulk = int(sys.argv[2])

legsize = 10

ion_list = ['HI', 'MgII', 'CIV', 'OVI']
galID_list = ['D9o2', 'D9q', 'D9m4a']
labels = ['dwSN', 'dwALL_1', 'dwALL_8']
mark=['v','o','x']
colors=['fuchsia', 'blue', 'black']
line=['-','.']


CL = 0.8413
tol = 1e-4
pmin = 0.0
pmax = 1.0
iter = 1000

# Define the EW cuts to use
ewcut = float(sys.argv[1])

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


# Location of master sysabs files
master_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/ALLfiles/masters/'

for ion in ion_list:
    i = 0
#    print ''
    for galID in galID_list:
        allfile = galID+'.'+ion+'.ALL.sysabs.large.master'
#        print ''
        absimpact=[]
        covering=[]
        imp=[]
        horizerrneg=[]
        horizerrpos=[]
        verterrneg=[]
        verterrpos=[]

        # Read in the allfile
        data = np.loadtxt(master_loc+allfile, skiprows=1)

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
            for k in range(0,len(data[:,0])):
                impact = data[k,22]
                absimpact.append(impact)
                width = data[k][5]
                
                if impact>minrad and impact<maxrad:
                    if width>ewcut:
                        hit += 1
                        xaxispoint += impact
                        total+=1
                        if minrad==0.0:
                            zerobinabs.append( data[k][0] )
                    else:
                        total+=1
            
            fraction = hit/total
            covering.append(fraction)
#            print ion+'\t'+galID+'\tRange: {0:0.1f} - {1:0.1f}\tHit: {2:d} \t Total: {3:d}\tCf: {4:.2f}'.format(minrad, maxrad, int(hit), int(total), fraction)

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

            
        if bulk==0:
            plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos],label=labels[i],marker=mark[i],c=colors[i],linestyle='none')
        else:
            subplotnum = 221+ion_list.index(ion)
            plt.subplot(subplotnum)
            plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos],label=labels[i],marker=mark[i],c=colors[i],linestyle='none')

        i+=1
        
    if bulk==0:
        plt.xlabel('D / R$_{vir}$')
        ylab = 'C$_f$ ['+ion+']'
        plt.ylabel(ylab)
        plt.ylim([-0.1,1.1])
        plt.xlim([min(Dmin),1.5])
        if ion=='HI':
            plt.legend(numpoints=1,loc='lower right',prop={'size':15}, frameon=False)
        else:
            plt.legend(numpoints=1,loc='upper right',prop={'size':15}, frameon=False)
        plt.gca().get_frame().set_linewidth(2) 
    
        plt.savefig('master_'+ion+'_{0:0.1f}mA_cover.eps'.format(ewcut*1000))
        plt.savefig('master_'+ion+'_{0:0.1f}mA_cover.pdf'.format(ewcut*1000))
        plt.savefig('master_'+ion+'_{0:0.1f}mA_cover.jpg'.format(ewcut*1000))
        plt.cla()
        plt.clf()

    else:
        subplotnum = 221+ion_list.index(ion)
        plt.subplot(subplotnum)
        plt.xlabel('D / R$_{vir}$')
        ylab = 'C$_f$ ['+ion+']'
        plt.ylabel(ylab)
        plt.ylim([-0.1,1.1])
        plt.xlim([min(Dmin),1.5])
        if ion=='HI':
            plt.legend(numpoints=1,loc='lower right',prop={'size':legsize}, frameon=False)
        else:
            plt.legend(numpoints=1,loc='upper right',prop={'size':legsize}, frameon=False)
        plt.gca().get_frame().set_linewidth(2) 

if bulk!=0:
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    s = 'master_{0:0.1f}mA_cover_bulk.eps'.format(ewcut*1000)
    print s
    
    plt.savefig(s, bbox_inches='tight')
