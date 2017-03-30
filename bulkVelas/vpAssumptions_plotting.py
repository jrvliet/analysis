
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white')


loc = '/home/jacob/research/code/analysis/bulkVelas/'
filename = 'vpAssumptions_vela.h5'
fname = loc+filename

df = pd.read_hdf(fname,'data')

galNums = range(21,30)
expns = np.arange(0.200,0.500,0.01)
expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]
ions = 'HI MgII CIV OVI'.split()

fields = 'loc.spread loc.sep length numCells'.split()
props = 'vlos t n Z'.split()
for prop in props:
    fields.append('{0:s}.mean'.format(prop))
    fields.append('{0:s}.std'.format(prop))



for field in fields:

    fig,axes = plt.subplots(2,2,figsize=(10,10))
    
    for ion,ax in zip(ions,axes.flatten()):

        for galNum in galNums:

            y = df[ion,field].ix[galNum]
            ax.plot(expns,y,label=galNum)

        ax.set_xlabel('Expn')
        ax.set_ylabel(field)
        ax.set_title(ion)
    fig.tight_layout()
    s = 'vpAssumptions_{0:s}.png'.format(field.replace('.','_'))
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)
            

