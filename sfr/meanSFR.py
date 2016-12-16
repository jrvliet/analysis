# coding: utf-8
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

loc = '/home/jacob/research/code/analysis/sfr/'
filename = 'vela2b-27_sfr.csv'
fname = loc+filename

df = pd.read_csv(fname)

window = 5
rollMeanSFR = df['sfr'].rolling(window=window,center=False).mean()

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(df['a'],rollMeanSFR)
ax.set_xlabel('Exapnsion Parameter')
ax.set_ylabel('Rolling Mean SFR')
ax.set_title('Window = {0;d}'.format(window))
fig.savefig('rollingMeanSFR.png',bbox_inches='tight',dpi=300)
plt.close(fig)


window = 5
rollMeanSFR = df['sfr'].rolling(window=3,center=False).mean()
fig,ax = plt.subplots(1,1,figsize=(5,5))
plt.plot(df['a'],rollMeanSFR)
plt.savefig('rollingMeanSFR.png')
plt.show()
df[df['sfr']<4]['a']
df['sfr'].mean()
df[df['sfr']<6]['sfr'].mean()
df[df['sfr']<3]['a']
df[df['sfr']<4]['a']
get_ipython().magic(u'save meanSFR.py')
get_ipython().magic(u'save meanSFR')
get_ipython().magic(u'save')
get_ipython().magic(u'save meanSFR 1-32')
