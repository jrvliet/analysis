
'''
Plots results of massContained.py
Only plot Rvir and proper kpc
'''


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

def mkLines(ax):
    ymin,ymax = ax.get_ylim()
    points = [0.29, 0.32, 0.35]
    points = [2.448,2.125,1.857,1.]
    for point in points:
        ax.vlines(point,ymin,ymax,linestyle='dashed',color='red')
    ax.set_ylim([ymin,ymax])

dataloc = '/home/jacob/research/code/analysis/inflows/'
filename = 'massContainedSimple.h5'

fractions = '25 50 75 90'.split()
masses = ['m{0:s}'.format(i) for i in fractions]
rvirs = ['rho{0:s}Rvir'.format(i) for i in fractions]
proper = ['rho{0:s}kpc'.format(i) for i in fractions]
comoving = ['rho{0:s}cokpc'.format(i) for i in fractions]

fields = [rvirs,proper,masses]


fname = dataloc+filename
store = pd.HDFStore(fname)

fig1,ax1 = plt.subplots(1,1,figsize=(5,5))
fig2,ax2 = plt.subplots(1,1,figsize=(5,5))
fig3,ax3 = plt.subplots(1,1,figsize=(5,5))
axes = [ax1,ax2,ax3]
figs = [fig1,fig2,fig3]
markers = 'o ^ s * x d'.split()
markers = 'o ^ s'.split()

labels = ['$-3.0<log(n_H)<-2.5$', '$-3.5<log(n_H)<-3.0$', 
          '$-4.0<log(n_H)<-3.5$', '$-4.0<log(n_H)<-4.5$', 
          '$-5.0<log(n_H)<-4.5$', '$-5.0<log(n_H)<-5.5$']
labels = ['$-3.5<log(n_H)<-2.5$',
          '$-4.5<log(n_H)<-3.5$',
          '$-5.5<log(n_H)<-4.5$']

# Loop over density cuts
for dfLabel,marker,label in zip(store,markers,labels):

    print(dfLabel)
    df = store[dfLabel]
    
    df['redshift'] = 1./df['a'] - 1

    for i,(field,ax) in enumerate(zip(fields,axes)):
        print('\t',field)
        ax.plot(df['redshift'],df[field[-1]],marker=marker,label=label)

for ax in axes:
    ax.invert_xaxis()
    ax.set_xlabel('Redshift')
    mkLines(ax)

ax1.legend(loc='lower left',fontsize='x-small')
ax2.legend(loc='upper left',fontsize='x-small')
ax3.legend(loc='center left',fontsize='x-small')

ax1.set_ylabel('$r_{90}$ [Rvir]')
ax2.set_ylabel('$r_{90}$ [pkpc]')
ax3.set_ylabel('Fraction of Mass Contained')

s1 = 'filamentRadius_rvir_simple.pdf'
s2 = 'filamentRadius_pkpc_simple.pdf'
s3 = 'filamentRadius_masses_simple.pdf'
fig1.savefig(s1,bbox_inches='tight',dpi=300)
fig2.savefig(s2,bbox_inches='tight',dpi=300)
fig3.savefig(s3,bbox_inches='tight',dpi=300)
plt.close(fig1)
plt.close(fig2)
plt.close(fig3)

store.close()

