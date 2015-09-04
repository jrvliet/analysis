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


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font',**font)

binnum = 20.0

# Filein name is structured as:
#   combined.(ion)[a,b].(expansion factor).sysabs

fileout = 'ew_v_D_combined_convert.pdf'
sym=['o','s']
pp = PdfPages(fileout)

numfiles = len(sys.argv)-1

for k in range(0,numfiles):

    filein = sys.argv[k+1]

    ionstring = filein.split('.')
    print ionstring
    ion = ionstring[2]
    print ion

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
    
    plt.plot(imp, meanew, 's', label=ion)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Impact Parameter (Rvir)')
ylab = "EW ($\AA$)"
plt.ylabel(ylab)
plt.ylim([0.2,3])
plt.legend(loc=1, frameon=False, numpoints=1)
plt.savefig(pp,format='pdf')
pp.close()
sys.exit()


hist, histedges = np.histogram(ew, bins=10, density=True)
print hist
print len(hist)
print len(histedges[1:])
print histedges
print histedges[1:]
#plt.plot(histedges[1:],hist,'sk')
plt.plot(impact, ew, marker='.', linestyle='none',label=gal)
#plt.hist(data[:,1].data[:,5],bins=20)
plt.xlabel('Impact Parameter (kpc)')
ylab = "EW ($\AA$)"
plt.ylabel(ylab)
plt.title(gal)
plt.yscale('log')
plt.xscale('log')
ymin, ymax = plt.ylim()
plt.xlim([1,100])
plt.savefig(pp,format='pdf')
plt.clf()
plt.cla()

pp.close()
