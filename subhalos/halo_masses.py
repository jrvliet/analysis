
'''
Plots evolution of mvir, stellar, gas mass of host halo
'''

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

dataloc = '/home/jacob/research/velas/vela2b/vela27/subhalos/'
filename = 'halos_{0:.3f}.txt'
expns = np.arange(0.200,0.550,0.01)

mvir, mgas, mstar, mdm = [],[],[],[]

for a in expns:

    fname = dataloc+filename.format(a)
    
    data = np.loadtxt(fname,skiprows=2)

    mvir.append(np.log10(data[0,7]))
    mdm.append(np.log10(data[0,11]))
    mgas.append(np.log10(data[0,12]))
    mstar.append(np.log10(data[0,13]))


fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))

ax1.plot(expns,mvir)
ax2.plot(expns,mdm)
ax3.plot(expns,mgas)
ax4.plot(expns,mstar)

for ax in [ax1,ax2,ax3,ax4]:
    ax.set_xlabel('Expansion Parameter')
    ymin,ymax = ax.get_ylim()
    ax.vlines(0.29,ymin,ymax,linestyle='dashed')
    ax.vlines(0.32,ymin,ymax,linestyle='dashed')
    ax.vlines(0.35,ymin,ymax,linestyle='dashed')


ax1.set_ylabel('Virial Mass')
ax2.set_ylabel('Dark Matter Mass')
ax3.set_ylabel('Gas Mass')
ax4.set_ylabel('Stellar Mass')




fig.tight_layout()
fig.savefig('vela2b-27_massEvolution.png',bbox_inches='tight',dpi=300)


