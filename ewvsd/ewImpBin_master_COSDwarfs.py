#!/usr/bin/python
#
# Filename: plotWvD_master_COSDwarfs.py
# Date: 04/09/14
#
# Usage:
#  python plotWvD_master_COSDwarfs.py 
#
# Plots the time evolution of a single ion
# Compares the different feedback perscriptions
# Bulk flag, if =1, will plot all four ions on one plot


import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv)!=3:
    bulk = 0
else:
    bulk = sys.argv[2]


galID_list = ['D9o2', 'D9q', 'D9m4a']
labels = ['dwSN', 'dwALL_1', 'dwALL_8']
ion_list = ['CIV']
mark=['v','o','x']
colors=['fuchsia', 'blue', 'black']
line=['-','.']

ewcut = 0.1
record = 1

# Define the radial bins to use
nbins = 15
binsize = 1.5/nbins

Dmax = []
Dmin = [0.0]
for j in range(0,nbins):
    Dmax.append( (j+1)*binsize )
    Dmin.append( (j+1)*binsize )
del Dmin[-1]
    
#for ion in ion_list:
ion = 'CIV'
outputs = []
for i in range(0,len(galID_list)):
    galID = galID_list[i]
    
    # Read in the master ALL file
    master_file = galID+'.'+ion+'.ALL.sysabs.large.master'
    data = np.loadtxt(master_file,skiprows=1)
    
    # Loop over the different impact parameter bins
    imp = []
    horizerrneg=[]
    horizerrpos=[]
    verterrneg=[]
    verterrpos=[]
    ewmean = []
    
    for j in range(0,len(Dmax)):
        minrad = Dmin[j]
        maxrad = Dmax[j]
        
        hit = 0.0
        xaxispoint = 0.0
        ewtotal = []
        
        # Loop over all the lines in the .sysabs file to get
        # the number of lines with impact parameters in this 
        # bin
        for k in range(0,len(data[:,0])):
            impact = data[k,22]
            ew = data[k,5]
            ew_err = data[k,6]
            if impact<maxrad and impact>minrad and ew_err!=-1.0 and ew>ewcut:
                hit += 1
                xaxispoint += impact
                ewtotal.append(ew)

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
        if standdev<=np.mean(ewtotal):
            verterrpos.append(standdev)
            verterrneg.append(standdev)
        else:
            verterrpos.append(standdev)
            verterrneg.append(np.mean(ewtotal)-1e-5)

    plt.errorbar(imp, ewmean, xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos], c=colors[i], marker=mark[i], label=labels[i], linestyle='none')


# Read in COS-Dwarfs Data if ion is CIV
cosdwarfs_file = '/home/matrix3/jrvander/data/obs_data/cos_dwarfs/cosdwarfs_detections.dat'
f  = open(cosdwarfs_file)
f.readline()
f.readline()
f.readline()

imp, rvir, ew, ew_err = [], [], [], []
for line in f:
    l = line.split()
    imp.append(float(l[7]))
    rvir.append(float(l[8]))
    
    ew.append(float(l[13])/1000.)
    ew_err.append(float(l[15])/1000.)
    
for i in range(0,len(imp)):
    imp[i] = imp[i] / rvir[i]

plt.errorbar(imp, ew, yerr=[ew_err, ew_err], markersize=5, marker='s', mec='black', c='cyan', ecolor='k', linestyle='none', label='COS-Dwarfs')




plt.xlabel('D / R$_{vir}$')
ylab = 'EW$_{CIV}$ [$\AA$]'
plt.ylabel(ylab)
plt.yscale('log')
plt.xlim([-0.05,1.5])
plt.ylim([5e-2, 2.0])
plt.legend(numpoints=1,loc=1,prop={'size':15}, frameon=False)
    
s = 'master_CIV_ewspace_COSDwarfs_cut_{0:d}mA.eps'.format(int(ewcut*1000))
print s

plt.savefig(s, bbox_inches='tight')
    
    
    


    
