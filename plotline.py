#!/usr/bin/python
#
# Filename: plotline.py
# 
# Plots the spectrum of lines created by specsynth
#
# Example Usage: 
#   python /home/matrix3/jrvander/code/plotline.py MW9a.MgII.los0001.MgII2796.spec MW9a.MgII.los0001.MgII2803.spec

import numpy as np
import pylab as plt
import sys

numfiles = len(sys.argv)-1
col = ['r', 'b']

peak = []
for i in range(0,numfiles):
    file = sys.argv[i+1]

    los = file.split('.')[2]
    line = file.split('.')[3]

    data = np.loadtxt(file)

    plt.step(data[:,1],data[:,2], label=line, color=col[i])

    peak.append(max(data[:,2])) 

plt.xlabel('Velocity (km/s)')
plt.ylabel('$I_{\lambda} / I^0_{\lambda}$')

if max(peak)<1.05:
    plt.ylim(ymax=1.05)

if 'HI' in file:
    xmin = -5000
    xmax = 5000
else:
    xmin = -1000
    xmax = 1000

plt.xlim([xmin, xmax])

plt.hlines([0,1], xmin, xmax, colors = 'k', linestyles = 'dashed')
plt.legend(loc=3, frameon=False)
plt.savefig('line_spectra.pdf')

