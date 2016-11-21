
'''
Plots results of massContained.py
'''


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

def mkLines(ax,ymin,ymax):
    points = [0.29, 0.32, 0.35]
    for point in points:
        ax.vlines(point,ymin,ymax,linestyle='dashed',color='red')

dataloc = '/home/jacob/research/code/analysis/inflows/'
filename = 'massContained.h5'

fractions = '25 50 75 90'.split()
masses = ['m{0:s}'.format(i) for i in fractions]
rvirs = ['rho{0:s}Rvir'.format(i) for i in fractions]
proper = ['rho{0:s}kpc'.format(i) for i in fractions]
comoving = ['rho{0:s}cokpc'.format(i) for i in fractions]

fields = [masses,rvirs,proper,comoving]


fname = dataloc+filename
store = pd.HDFStore(fname)

fig2,axes2 = plt.subplots(2,2,figsize=(10,10))
properLims = [0,250]
comoveLims = [0,700]

# Loop over density cuts
for dfLabel in store:

    df = store[dfLabel]
    
    fig,axes = plt.subplots(2,2,figsize=(10,10))
    # Loop over each plot, which represents a different quantity such as mass,
    # rho in Rvir, pkpc, or ckpc
    for i,(field,ax) in enumerate(zip(fields,axes.flatten())):
        # Loop over each fraction cut
        for f,fraction in zip(field,fractions):
            ax.plot(df['a'],df[f],label=fraction)

    # Apply labels to plots and add vertical lines for each major event
    ylabs = ['Mass', r'$\rho$ [Rvir]', r'$\rho$ [pkpc]', r'$\rho$ [ckpc]']
    for i,ax in enumerate(axes.flatten()):
        ax.set_xlabel('a')
        ax.set_ylabel(ylabs[i])
    axes[1,0].set_ylim(properLims)
    axes[1,1].set_ylim(comoveLims)
    for ax in axes.flatten():
        ymin,ymax = ax.get_ylim()
        mkLines(ax,ymin,ymax)
    
    ax.legend()
    fig.tight_layout()
    s ='massContained_{0:s}.png'.format(dfLabel[1:].replace('(',
                                    '').replace(')',''))
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


    # Make a plot similar to above, but with only the 90% cut in it to directly
    # compare the affect of density cuts
    for i,(field,ax) in enumerate(zip(fields,axes2.flatten())):
        ax.plot(df['a'],df[field[-1]],label=dfLabel[1:])


for i,ax in enumerate(axes2.flatten()):
    ax.set_xlabel('a')
    ax.set_ylabel(ylabs[i])
    ymin,ymax = ax.get_ylim()
    mkLines(ax,ymin,ymax)
ax.legend()
fig2.tight_layout()
s = 'massContained_denseCuts.png'
fig2.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig2)

store.close()

