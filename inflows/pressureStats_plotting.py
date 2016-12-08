


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



dataloc = '/home/jacob/research/code/analysis/inflows/'
filename = 'pressureStats_noCGM.h5'

fname = dataloc+filename

df = pd.read_hdf(fname,'data')
df['a'] = np.arange(0.200,0.550,0.01)


fig,ax = plt.subplots(1,1,figsize=(5,5))
for name in df.columns:
    if 'Mean' in name:
        ax.plot(df['a'],df[name],label=name.strip('Mean')+'Z')

ax.set_xlabel('Expansion Parameter')
ax.set_ylabel('Mean Pressure')
ax.set_yscale('log')
ax.legend(loc='best')
ymin,ymax = ax.get_ylim()
ax.vlines(0.29,ymin,ymax,linestyle='dashed')
ax.vlines(0.32,ymin,ymax,linestyle='dashed')
ax.vlines(0.35,ymin,ymax,linestyle='dashed')
s = 'pressureStats_mean.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)

fig,ax = plt.subplots(1,1,figsize=(5,5))
for name in df.columns:
    if 'Std' in name:
        ax.plot(df['a'],df[name],label=name.strip('Std')+'Z')

ax.set_xlabel('Expansion Parameter')
ax.set_ylabel('Std Pressure')
ax.set_yscale('log')
ymin,ymax = ax.get_ylim()
ax.vlines(0.29,ymin,ymax,linestyle='dashed')
ax.vlines(0.32,ymin,ymax,linestyle='dashed')
ax.vlines(0.35,ymin,ymax,linestyle='dashed')
ax.legend(loc='best')
s = 'pressureStats_std.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)



