
'''
Plot a histgram of the standard deviations reported by clumping.py
'''

import numpy as np
import matplotlib.pyplot as plt
import json

filename = 'stddev.dat'
numbins = 50

# Read in header
#with open(filename, 'r') as f:
#    header = f.readline().strip()

with open(filename, 'r') as f:
    data = json.load(f)

# Parse the ion names
#ions = [ s for s in  header.split()]

# Read in the rest of the data
#data = np.loadtxt(filename, skiprows=1)

f = open('stddev.out', 'w')
header = '{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\n'.format(
    'ion', 'Min', 'Max', 'Mean', 'Std Dev')
f.write(header)
formatString = '{0:s}\t{1:.3f}\t{2:.3f}\t{3:.3f}\t{4:.3f}\n'

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

#ions = ['HI', 'MgII', 'CIV', 'OVI']

ions = data[0]
for i in range(0,len(ions)):

    ion = ions[i]
    d = data[i+1]
    ax1.hist(d, bins=numbins, range=(0,1), histtype='step', label=ion)
    ax2.hist(d, bins=numbins, range=(0,1), histtype='step', log=True, label=ion)

    f.write(formatString.format(ion, min(d), max(d), np.mean(d), np.std(d)))

ax1.legend(frameon=False)
ax1.set_xlabel('Std Dev of Normalized Cell Loc')
ax1.set_ylabel('Counts')
fig1.savefig('stdDevHist.png', bbox_inches='tight')

ax2.legend(frameon=False)
ax2.set_xlabel('Std Dev of Normalized Cell Loc')
ax2.set_ylabel('Log Counts')
fig2.savefig('stdDevHist_log.png', bbox_inches='tight')

f.close()


    


