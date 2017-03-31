
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white')
import sys


loc = '/home/jacob/research/code/analysis/bulkVelas/'
filename = 'vpAssumptions_vela.h5'
fname = loc+filename

df = pd.read_hdf(fname,'data')

galNums = range(21,30)
galNums.remove(27)
expns = np.arange(0.200,0.550,0.01)
expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]
ions = 'HI MgII CIV OVI'.split()

fields = 'loc.spread loc.sep length numCells'.split()
props = 'vlos t n Z'.split()
for prop in props:
    fields.append('{0:s}.mean'.format(prop))
    fields.append('{0:s}.std'.format(prop))

logFields = 't n Z'.split()

loLimits = [0, 0, 0, 0, -250, 0, 4, 3.4, -3.5, -4, -5.6, -6]
hiLimits = [1, 0.25, 1400, 250, 250, 110, 6.5, 6.4, 2, 2, -4, -3.8]

for field,lo,hi in zip(fields,loLimits,hiLimits):
    print(field)

    fig,axes = plt.subplots(2,2,figsize=(10,10))
    
    for ion,ax in zip(ions,axes.flatten()):
        for galNum in galNums:
            x = expns
            y = df[ion,field].ix[galNum].astype(float)
            if field[0] in logFields:
                y = np.log10(y)
            #print(len(x),len(y))
            #print(x)
            #print(y)
            ax.plot(x,y,label=galNum)

        ax.set_ylim([lo,hi])
        ax.set_xlabel('Expn')
        ax.set_ylabel(field)
        ax.set_title(ion)
    axes[0,0].legend(loc='best')
    fig.tight_layout()
    s = 'vpAssumptions_{0:s}.png'.format(field.replace('.','_'))
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)
            

