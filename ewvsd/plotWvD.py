#!/usr/bin/python
#
# Filename: plotWvD.py
# Date: 17/12/13
#
# Usage:
#  python plotWvD.py <gal_props_file> <ALL file> 


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







import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages


import matplotlib


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font',**font)


# First sysarg is the name of the control file, usually gal_props.dat
controlfile = sys.argv[1]
print controlfile
gasfile, galID, rootname, expn, summaryLoc = read_control_file(controlfile)
mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z = read_summary(galID, expn, summaryLoc)


ion = sys.argv[1].split('.')[1]
#filein = ['MW9a.CIV.ALL.sysabs']

#fil = sys.argv[1]
#filein = [fil]

fileout = galID+'_{0:1.3f}_ewspace.pdf'.format(expn)
sym=['o','s']

pp = PdfPages(fileout)
i=0

data = np.loadtxt(sys.argv[2],skiprows=1)

#plt.plot(data[:,1],data[:,5],marker='.', linestyle='none',label=galID)

plt.xlabel('Impact Parameter (kpc)')
ylab = "EW ($\AA$)"
plt.ylabel(ylab)
plt.title(galID+' '+ion)
plt.yscale('log')
plt.xscale('log')
ymin, ymax = plt.ylim()

   
#    if i%2==1:
plt.savefig(pp,format='pdf')
plt.clf()
plt.cla()

i+=1

pp.close()
