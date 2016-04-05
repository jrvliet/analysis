
'''
Bins up the LOS velocity and calculates the standard deviation of all 
cell's position along that LOS that have that LOS velocity
'''

import numpy as np
import matplotlib.pyplot as plt
import json

ions = ['HI', 'MgII', 'CIV', 'OVI']

numbins = 50
vmax = 100.0
vmin = -100.0

binsize = (vmax-vmin)/numbins
vlo = []
vhi = []
for i in range(numbins):
    vlo.append(i*binsize + vmin)
    
print vlo

for ion in ions:

    filename = '{0:s}_vlos.dat'.format(ion)
    l, s, v = np.loadtxt(filename, skiprows=1, usecols=(0,1,2), unpack=True)
    
    stdev, xloc = [], []

    for i in range(0,numbins):
        loc = []
        lower = vlo[i]
        if i==numbins-1:
            upper = vmax
        else:
            upper = vlo[i+1]
        binmid = (upper+lower)/2.0
        for j in range(0,len(v)):
            if v[j]>lower and v[j]<=upper:
                loc.append(s[j]*l[j])

        stdev.append(np.std(loc))
        xloc.append(binmid)

    plt.plot(xloc, stdev, 'x-', label=ion)


plt.xlabel('LOS velocity')
plt.ylabel('Std Dev of LOS Position')
plt.legend(frameon=False)

plt.savefig('binnedVlos.png', bbox_inches='tight')


























        

