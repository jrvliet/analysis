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
import sys

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
    vmin,vmax = 0,4.3
    hnorm,xedges,yedges,binnumber = st.binned_statistic_2d(n,t,z,
                                statistic=stat,range=binrange,
                                bins=numbins)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h) 
    h = np.log10(h)
    mesh = ax.pcolormesh(xedges,yedges,h,vmin=vmin,vmax=vmax)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    
    return hnorm,xedges,yedges
    
kpc2km = 3.086e16
pc2cm = 3.086e18
mH = 1.6737e-24
mSun = 1.989e33
s2yr = 3.154e7

vmin,vmax = 0,4.3

rootloc = '/home/jacob/research/velas/vela2b/'
rootloc = '/mnt/cluster/abs/cgm/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i{2:d}/{3:s}/'
filename = '{0:s}.{1:s}.a{2:.3f}.i{3:d}.ALL.sysabs.h5'
filename = '{0:s}.{1:.3f}.{2:s}.i90.abs_cells.h5'

galNums = range(21,30)

ions = 'HI MgII CIV OVI'.split()

finalExpn = [0.530, 0.490, 0.490, 0.510, 0.460, 
             0.500, 0.500, 0.500, 0.510, 0.500]

expns = np.arange(0.200,max(finalExpn),0.01)
expns = np.arange(0.200,0.500,0.01)
expnLabels = ['a{0:d}'.format(int(a*1000)) for a in expns]

lo,hi = [],[]

count = 0
grandtotals = [ [] for i in range(len(ions))]

errlogFile = 'phase.errlog'
errlog = open(errlogFile,'w')
header = 'galID\texpn\tion\n'
errlog.write(header)
errs = '{0:s}\t{1:.3f}\t{2:s}\n'
for a,aLabel in zip(expns,expnLabels):
    #fig,axes = plt.subplots(len(galNums),len(ions),
    #                        figsize=(20,45),
    #                        sharex=True,
    #                        sharey=True)

    print(aLabel)
    #totals = [None]*len(ions)
    for i,(galNum,finala) in enumerate(zip(galNums,finalExpn)):
        loc = rootloc+'vela{0:d}/a{1:.3f}/i90/'.format(galNum,a)

        print('\t',galNum)
        try:
            galID,expn,redshift,mvir,rvir,inc = galaxyProps(loc)

            for j,ion in enumerate(ions):

                #ax = axes[i,j]
                
                cells = rootloc+subloc.format(galNum,a,inc,
                            ion)+filename.format(galID,a,ion)
                df = pd.read_hdf(cells,'data')
                try:
                    df['mass'] = (10**df['nH'])*mH*(df['cell_size']*pc2cm)**3 / mSun
 
                    numbins = 500
                    binrange = [[-8,3],[2,8]]
                    vmin,vmax = 0,4.3
                    stat = 'sum'
                    n = df['nH']
                    t = df['temperature']
                    z = df['mass']
                    h,xedges,yedges,binnumber = st.binned_statistic_2d(n,t,z,
                                                statistic=stat,range=binrange,
                                                bins=numbins)                   
                    h = np.rot90(h)
                    h = np.flipud(h)
                    #h,xedges,yedges = mkHist(ax,df['nH'],df['temperature'],
                    #                df['mass'],'count')
                    savename = '{0:s}_{1:s}_{2:d}_phase.dat'.format(ion,aLabel,galNum)
                    np.savetxt(savename,h)
#
#                    if i==0:
#                        totals[j] = h
#                    else:
#                        #print(type(totals[j],type(h))
#                        totals[j] = totals[j]+h
#
#                    lo.append(h.min())
#                    hi.append(h.max())
#
#                    grandtotals[j].append(h)
#                    count+=1
                except KeyError:
                    print(galID,a,ion)
                    errlog.write(errs.format(galID,a,ion))
                    continue

        except IOError:
            print(cells)
            continue

#    for ax,ion in zip(axes[0,:],ions):
#        ax.set_title(ion)
#    for ax,galNum in zip(axes[:,0],galNums):
#        ax.set_ylabel('vela2b-{0:d}\n$\log(T)$ [K]'.format(galNum))
#    fig.subplots_adjust(wspace=0,hspace=0)
#    s = 'vela2b_{0:s}_phase.png'.format(aLabel)
#    fig.savefig(s,bbox_inches='tight',dpi=300)
#    plt.close(fig)
#
#    # Plot Totals
#    fig,axes = plt.subplots(2,2,figsize=(10,10),sharex=True,sharey=True)
#    for ion,ax,total in zip(ions,axes.flatten(),totals):
#        ax.pcolormesh(xedges,yedges,total,vmin=vmin,vmax=vmax)
#        #mkHist(ax,total['nH'],total['temperature'],total['mass'],'mean')
#        binrange = [[-10,2],[2,8]]
#        ax.set_xlim(binrange[0])
#        ax.set_ylim(binrange[1])
#        ax.annotate(ion,(-8,3),xycoords='data')
#    fig.subplots_adjust(wspace=0,hspace=0)
#    fig.savefig('vela2b_{0:s}mean_phase.png'.format(aLabel),bbox_inches='tight',dpi=300)
#    plt.close(fig)

# Save totals to file
for ion,grand in zip(ions,grandtotals):
    for d,lab in zip(grand,expnLabels):
        filename = '{0:s}_{1:s}_phase.dat'.format(ion,aLabel)
        np.savetxt(filename,h)
sys.exit()

    
# Plot grand total and std
fig,axes = plt.subplots(2,2,figsize=(10,10),sharex=True,sharey=True)
cbarax = fig.add_axes([0.91, 0.125, 0.01, 0.775])
for ion,ax,grand in zip(ions,axes.flatten(),grandtotals):
    print(type(grand))
    grand = [10**i for i in grand]
    ave = np.average(grand,axis=0)
    ave = np.ma.masked_where(ave==0,ave)
    ave = np.log10(ave)
    mesh = ax.pcolormesh(xedges,yedges,ave,vmin=vmin,vmax=vmax)
    binrange = [[-10,2],[2,8]]
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.annotate(ion,(-8,3),xycoords='data')
    fig.subplots_adjust(wspace=0,hspace=0)
cbar = fig.colorbar(mesh,ax=cbarax,use_gridspec=True)
for ax in axes[-1,:]:
    ax.set_xlabel('$\log ( n_H $ [cm$^{-3}$] )')
for ax in axes[:,0]:
    ax.set_ylabel('$\log ( T $[K] ')
s = 'vela2b_phase_fullMean.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)

fig,axes = plt.subplots(2,2,figsize=(10,10),sharex=True,sharey=True)
cbarax = fig.add_axes([0.91, 0.125, 0.01, 0.775])
for ion,ax,grand in zip(ions,axes.flatten(),grandtotals):
    print(type(grand))
    grand = [10**i for i in grand]
    ave = np.std(grand,axis=0)
    ave = np.ma.masked_where(ave==0,ave)
    ave = np.log10(ave)
    mesh = ax.pcolormesh(xedges,yedges,ave,vmin=vmin,vmax=vmax)
    binrange = [[-10,2],[2,8]]
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.annotate(ion,(-8,3),xycoords='data')
    fig.subplots_adjust(wspace=0,hspace=0)
cbar = fig.colorbar(mesh,ax=cbarax,use_gridspec=True)
for ax in axes[-1,:]:
    ax.set_xlabel('$\log ( n_H $ [cm$^{-3}$] )')
for ax in axes[:,0]:
    ax.set_ylabel('$\log ( T $[K] ')
s = 'vela2b_phase_fullStd.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)

     

errlog.close()






print(min(lo),max(hi))

