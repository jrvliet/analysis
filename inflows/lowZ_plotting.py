
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


dataloc = '/home/jacob/research/code/analysis/inflows/'
filename = 'lowZProps_filament.h5'

fname = dataloc+filename
df = pd.read_hdf(fname,'data')

fields = df.columns[1:]

for field in fields:

    fig,ax = plt.subplots(1,1,figsize=(5,5))
    ax.plot(df['a'],df[field])
    ax.set_xlabel('Expn')
    ax.set_ylabel(field)
    
    s = 'lowZProps_{0:s}.png'.format(field)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)







