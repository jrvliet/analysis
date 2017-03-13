#!/usr/bin/python

'''
Bins up the phase space of 
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
pd.options.mode.chained_assignment = None

def galaxyProps(location):

    fname = location+'galaxy.props'
    with open(fname) as f:
        galID = f.readline().split()[1]
        expn = f.readline().split()[1]
        redshift = f.readline().split()[1]
        mvir = f.readline().split()[1]
        rvir = float(f.readline().split()[1])
        inc = int(float(f.readline().strip().split()[1]))

    return galID,expn,redshift,mvir,rvir,inc



def mkHist(ax,n,t,z,stat):

    numbins = 50
    binrange = [[-10,2],[2,8]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(n,t,z,
                                statistic=stat,range=binrange,
                                bins=numbins)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    
    h = np.log10(h)
    mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    
    return h.min(),h.max()
    
kpc2km = 3.086e16
pc2cm = 3.086e18
mH = 1.6737e-24
mSun = 1.989e33
s2yr = 3.154e7


rootloc = '/mnt/cluster/abs/cgm/vela2b/'
rootloc = '/home/jacob/research/velas/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i{2:d}/{3:s}/'
filename = '{0:s}.{1:s}.a{2:.3f}.i{3:d}.ALL.sysabs.h5'
filename = '{0:s}.{1:.3f}.{2:s}.i90.abs_cells.h5'

galNums = range(21,30)

ions = 'HI MgII CIV OVI'.split()

finalExpn = [0.530, 0.490, 0.490, 0.510, 0.460, 
             0.500, 0.500, 0.500, 0.510, 0.500]

expns = np.arange(0.200,max(finalExpn),0.01)
expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]

lo,hi = [],[]
for a,aLabel in zip(expns,expnLabels):
    fig,axes = plt.subplots(len(galNums),len(ions),
                            figsize=(20,45),
                            sharex=True,
                            sharey=True)

    for i,(galNum,finala) in enumerate(zip(galNums,finalExpn)):
        loc = rootloc+'vela{0:d}/a{1:.3f}/i90/'.format(galNum,a)

        try:
            galID,expn,redshift,mvir,rvir,inc = galaxyProps(loc)

            for j,ion in enumerate(ions):

                ax = axes[i,j]
                
                cells = rootloc+subloc.format(galNum,a,inc,
                            ion)+filename.format(galID,a,ion)
                df = pd.read_hdf(cells,'data')
                df['mass'] = (10**df['nH'])*mH*(df['cell_size']*pc2cm)**3 / mSun

                clo,chi = mkHist(ax,df['nH'],df['temperature'],
                                df['mass'],'count')
                lo.append(clo)
                hi.append(chi)

        except IOError:
            continue

    s = 'vela2b_a{0:s}_phase.png'.format(aLabel)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)
    break

print(min(lo),max(hi))

