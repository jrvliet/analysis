
# coding: utf-8

# In[14]:

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



loc = '/home/jacob/research/code/analysis/bulkVelas/'
filename = 'covering.h5'
fname = loc+filename

df = pd.read_hdf(fname,'data')


sns.set_style('white')
ions = 'HI MgII CIV OVI'.split()
for ion in ions:
    hi = df.xs(ion,level=1,axis=1).transpose()
    print(ion,hi.values.min(),hi.values.max())
    hi = hi.mask(hi==0)
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    sns.heatmap(hi,ax=ax,vmin=0,vmax=1,cmap='viridis',annot=False,fmt='0.1f')
    ax.set_title(ion)
    s = '{0:s}/covering_{1:s}.png'.format(loc,ion)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)




