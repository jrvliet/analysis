#!/usr/bin/python
#
# Filename: plotWvD.py
# Date: 17/12/13
#
# Usage:
#  python plotWvD.py <ion>
#
# Plots the time evolution of a single ion


def read_control_file(filename):

    f = open(filename)
    
    for i in range(0,5):
        line = f.readline()

    gasfile = line.split()[0]
    gasfile = gasfile.replace('"', '')
    galID = gasfile.split('_')[0]
    
    line = f.readline()
    rootname = line.split()[0]
    
    line = f.readline()
    aexpn = float(line.split()[0])
    
    for i in range(0,3):
        line = f.readline()
    summaryLoc = (line.split()[0]).replace('"','')

    f.close()

    return gasfile, galID, rootname, aexpn, summaryLoc

def read_summary(galID, aexpn, summaryLoc):

    import sys

    aexpn = float(aexpn)

    # Location of summaries
    f = open(summaryLoc + galID + '.dat')

    f.readline()
    f.readline()

    found = 0

    for line in f:
#    while(found==0):

        l = line.split()
        a = float(l[0])

        if a==aexpn:
            
            found = 1

            # Read in data
            mvir = float(l[2])
            rvir = float(l[3])
            a11 = float(l[4])
            a12 = float(l[5])
            a13 = float(l[6])
            a21 = float(l[7])
            a22 = float(l[8])
            a23 = float(l[9])
            a31 = float(l[10])
            a32 = float(l[11])
            a33 = float(l[12])
            vpec_x = float(l[13])
            vpec_y = float(l[14])
            vpec_z = float(l[15])
            
    f.close()
    if found==1:
        return mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z
    else:
        print 'ERROR: Could not find the expansion factor {0} in file: \n\t'.format(aexpn) + summaryLoc
        sys.exit()

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
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
import scipy.special as sc
import scipy.optimize as opt
import math
import matplotlib
import os


font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

#matplotlib.rc('font',**font)


expn_list = ['0.625', '0.650', '0.675', '0.700', '0.725']
ion = sys.argv[1]
mark=['v','o','*','s','x']
colors=['red', 'darkorange', 'green', 'cyan', 'blue']
line=['-','.']


# Define the radial bins to use
nbins = 15
binsize = 1.5/nbins

# Define an EW cut
ewcut = 50.0e-3   # Estimated from COS-Dwarfs presentations


# Loop over the different expansion factors
for i in range(0,len(expn_list)):
    
    # Move into that directory
    os.chdir('./a'+expn_list[i]+'/'+ion+'/')


    # Read in galaxy properties
    controlfile = 'gal_props_'+ion+'.dat'
    gasfile, galID, rootname, expn, summaryLoc = read_control_file(controlfile)
    mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z = read_summary(galID, expn, summaryLoc)


    Dmax = []
    Dmin = [0.0]
    for j in range(0,nbins):
        Dmax.append( (j+1)*binsize )
        Dmin.append( (j+1)*binsize )
    del Dmin[-1]


    fileout = galID+'_{0:1.3f}_ewspace.pdf'.format(expn)
    sym=['o','s']

    # Read in the ALL file
    allfile = galID+'.'+ion+'.ALL.sysabs'
    data = np.loadtxt(allfile,skiprows=1)
#lines = np.loadtxt('lines.info',skiprows=2)

    # Loop over the different impact parameter bins
    imp = []
    horizerrneg=[]
    horizerrpos=[]
    verterrneg=[]
    verterrpos=[]
    ewmean = []
    for j in range(0,len(Dmax)):
        minrad = Dmin[j]*rvir
        maxrad = Dmax[j]*rvir
        
        #    print 'Minrad: ', minrad
        #    print 'Maxrad: ', maxrad
        hit = 0.0
        xaxispoint = 0.0
        ewtotal = []

        # Loop over all the lines in the .sysabs file to get
        # the number of lines with impact parameters in this 
        # bin
        for k in range(0,len(data[:,0])):
            impact = data[k,1]
            ew = data[k,5]
            ew_err = data[k,6]
            if impact<maxrad and impact>minrad and ew_err!=-1.0 and ew>ewcut:
                hit += 1
                xaxispoint += impact
                ewtotal.append(ew)
#    print ewtotal
#    print 'Hit: ',hit
#    print ''
        # Take the mean of ewtotal
        ewmean.append(np.mean(ewtotal))
    

        # Determine the x location of the point. If hit>0, then the location is
        # the average impact parameter in that bin. If hit=0, then the location
        # is the midpoint of the bin
        if hit > 0.0:
            imp.append(xaxispoint/hit)  # location of point on x-axis
        else:
            imp.append((maxrad-minrad)/2.0 + minrad)
    
        # Determine error bars
        #   Horizontal error bars are width of bin
        horizerrneg.append(imp[j]-minrad)
        horizerrpos.append(maxrad-imp[j])
           
        #   Vertical error bars are found using the standard deviation
        standdev = np.std(ewtotal)
        verterrpos.append(standdev)
        verterrneg.append(standdev)

    imp = [k/rvir for k in imp]
    horizerrneg = [k/rvir for k in horizerrneg]
    horizerrpos = [k/rvir for k in horizerrpos]
    z = math.fabs(1.0/float(expn_list[i]) - 1.0)
    plt.errorbar(imp, ewmean, xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos], c=colors[i], marker=mark[i], label='z=:{0:1.2f}'.format(z), linestyle='none')

    os.chdir('../..')


#plt.xlabel('Impact Parameter (Rvir)')
plt.xlabel('$\\rho$ / R$_{vir}$')
ylab = 'EW$_{'+ion+'}$ [$\AA$]'
plt.ylabel(ylab)
#plt.title(galID +'; '+ion)
#plt.title(ion)
plt.yscale('log')
ymin, ymax = plt.ylim()
plt.xlim([-0.05,1.5])
if ion=='HI':
    plt.ylim([0.01, 10.0])
else:
    plt.ylim([0.001, 1.0])
plt.legend(numpoints=1,loc=1,prop={'size':15}, frameon=False)

s = galID+'_'+ion+'_ewspace_cut_{0:d}mA.pdf'.format(int(ewcut*1000))
print s

plt.savefig(s)
