

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('vela2b-27_sfr.csv')
df['redshift'] = 1./df['a'] - 1

fig,ax1 = plt.subplots(1,1,figsize=(5,5))
ax1.invert_xaxis()
ax2 = ax1.twinx()

ln1 = ax1.plot(df['redshift'],df['sfr'],label='SFR',color='blue')
ln2 = ax2.plot(df['redshift'],df['ssfr'],label='sSFR',color='red')
ax1.set_xlabel('Redshift')
ax1.set_ylabel('SFR')
ax2.set_ylabel('sSFR')

lns = ln1+ln2
labs = [l.get_label() for l in lns]
ax1.legend(lns,labs,loc='best')

s = 'sfrPlay.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)



quants = np.arange(0.2,1.0,0.1)
fig,axes = plt.subplots(2,4,figsize=(20,10))
for quant,ax in zip(quants,axes.flatten()):
    sfrCut = df['sfr'].quantile(quant)
    a = df['a'][df['sfr']>sfrCut]
    sfr = df['sfr'][df['sfr']>sfrCut]
    
    ax.plot(df['a'],df['sfr'],color='blue')
    ax.plot(a,sfr,color='red',marker='s',linestyle='none')
    
    ax.set_xlabel('Expansion Parameter')
    ax.set_ylabel('SFR')
    ax.set_title('Cut = {0:.1f}'.format(quant))

fig.tight_layout()
s = 'sfrQuantiles.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)



sigma = df['sfr'].std()
mean = df['sfr'].mean()
median = df['sfr'].median()
fig,axes = plt.subplots(1,4,figsize=(20,5))
multiples = np.arange(0.5, 2.5, 0.5)
print(sigma,mean,median)
df['sfrMean'] = df['sfr']-median
for mult,ax in zip(multiples,axes.flatten()):

    sfrMult = mult*sigma
    cutInd = (df['sfr']>mean+sfrMult)
    a = df['a'][cutInd]
    sfr = df['sfr'][cutInd]

    print('\t',sfrMult,len(a))
    ax.plot(df['a'],df['sfr'],color='blue')
    ax.plot(a,sfr,color='red',marker='s',linestyle='none')
    
    ax.set_xlabel('Expansion Parameter')
    ax.set_ylabel('SFR')
    ax.set_title('Multiple = {0:.1f}$\sigma$'.format(mult))
    
fig.tight_layout()
s = 'sfrSigma.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)

















