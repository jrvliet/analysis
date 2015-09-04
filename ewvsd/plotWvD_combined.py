#!/usr/bin/python
#
# Filename: plotWvD_combined.py
# Date: 12/31/12
#
# Description:
#  Read in results from Chris' bindata-logfreq program and plot it
#  Works for the combined results from several galaxies


import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
from pylab import *

def log_10_product(x, pos):
    """The two args are the value and tick position.
    Label ticks with the product of the exponentiation"""
    if x>=1:
        str = '%li' % (x)
    else:
        str = '%0.1f' % (x)
    return str

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font',**font)

binnum = 20.0

# Filein name is structured as:
#   combined.[a,b].ion.(expansion factor).sysabs

fileout = 'ew_v_D_combined.pdf'
sym=['o','s']
col=['k','r']


numfiles = len(sys.argv)-1

sf = sys.argv[1].split('.')[1]
ion = sys.argv[1].split('.')[2]
expansion = sys.argv[1].split('.')[3]+'.'+sys.argv[1].split('.')[4]

fileout = 'ew_v_d_combined_'+ion+'_'+expansion+'.pdf'
print 'Output file: ',fileout
pp = PdfPages(fileout)

ax = subplot(111)
for k in range(0,numfiles):

    filein = sys.argv[k+1]
    ionstring = filein.split('.')
    ion = ionstring[2]
    sf = ionstring[1]

    data = np.loadtxt(filein,skiprows=1)

    impact = data[:,1]
    ew = data[:,5]
    maxb = max(impact)
    minb = min(impact)
    step = (maxb-minb)/binnum
    
    meanew = []
    imp = []
    for i in range(0,int(binnum)):
        low = step*i + minb
        high = low + step
        imp.append((low+high)/2.0)
        ewsum = 0.0
        count = 0.0
        for j in range(0,len(impact)):
            if impact[j]>low and impact[j]<high:
                ewsum += ew[j]
                count += 1
        meanew.append(ewsum / count)
 
    plt.plot(imp, meanew, 's', label=sf)

plt.xscale('log')
plt.yscale('log')

formatter = FuncFormatter(log_10_product)

ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)
plt.grid(b=True, which='major', axis='both', ls='--')
plt.tick_params(axis='both', which='minor', length=5, width=1)
plt.xlabel('Impact Parameter b (kpc)')
ylab = "Rest Equivalent Width $W_0$ ($\AA$)"
plt.ylabel(ylab)
plt.ylim([0.02,10])
plt.xlim([1,250.0])
plt.legend(loc=1, frameon=False, numpoints=1)
plt.savefig(pp,format='pdf')
pp.close()
sys.exit()
