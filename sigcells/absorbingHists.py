
'''
Plot a histgram of the standard deviations reported by clumping.py
'''

import numpy as np
import matplotlib.pyplot as plt
import glob
import json
import sys

# Get the filenames available
files = glob.glob('std*.dat')

numbins = 50

formatString = '{0:s}\t{1:.3f}\t{2:.3f}\t{3:.3f}\t{4:.3f}\n'
header = '{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\n'.format(
    'ion', 'Min', 'Max', 'Mean', 'Std Dev')

fig1, ((ax11,ax12),(ax13,ax14)) = plt.subplots(2,2,figsize=(10.2,10.2))
fig2, ((ax21,ax22),(ax23,ax24)) = plt.subplots(2,2,figsize=(10.2,10.2))

axes1 = [ax11,ax12,ax13,ax14]
axes2 = [ax21,ax22,ax23,ax24]

for j, filename in enumerate(files):
    
    prop = filename.split('std')[1].split('.')[0]
    print filename, prop

    ax1 = axes1[j]
    ax2 = axes2[j]

    with open(filename, 'r') as f:
        data = json.load(f)

    with open(filename.replace('dat','out'), 'w') as f:
        f.write(header)

        ions = data[0]
        for i in range(0,len(ions)):

            
            ion = ions[i]
            d = data[i+1] 
            if prop=='Temperatures':
                d = np.log10(d)
                dmin = 0
                dmax = 2
            elif prop=='Densities':
                d = np.log10(d)
                dmin = 0
                dmax = 2
            elif prop=='Metals':
                d = np.log10(d)
                dmin = 0
                dmax = 1 
            elif prop=='Location':
                dmin = 0
                dmax = 150
            else:
                print 'Unknown property ',prop
                sys.exit()
                   
#            ax1.hist(d, bins=numbins, histtype='step', range=[dmin,dmax], label=ion)
#            ax2.hist(d, bins=numbins, histtype='step', range=[dmin,dmax], log=True, label=ion)
            ax1.hist(d, bins=numbins, histtype='step', label=ion)
            ax2.hist(d, bins=numbins, histtype='step', log=True, label=ion)

            f.write(formatString.format(ion, min(d), max(d), np.mean(d), np.std(d)))

        ax1.legend(frameon=False)
        ax1.set_xlabel('Std Dev of {0:s}'.format(prop))
        ax1.set_ylabel('Counts')

        ax2.legend(frameon=False)
        ax2.set_xlabel('Std Dev of {0:s}'.format(prop))
        ax2.set_ylabel('Log Counts')


fig1.tight_layout()
fig2.tight_layout()
    
plotname = 'stdDevHist_quad.pdf'
fig1.savefig(plotname, bbox_inches='tight')
plotname = plotname.replace('.pdf','_log.pdf')
fig2.savefig(plotname, bbox_inches='tight')


