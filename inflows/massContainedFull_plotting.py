
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def mkLines(ax):
    ymin,ymax = ax.get_ylim()
    points = [2.57,2.13,1.86]
    for p in points:
        ax.vlines(p,ymin,ymax,color='black',linestyle='--')
    ax.set_ylim([ymin,ymax])

loc = '/home/jacob/research/code/analysis/inflows/'
filename = 'filamentMassSize_simple.h5'
fname = loc+filename

df = pd.read_hdf(fname,'data')

denseLabels = 'rare mid dense'.split()
tempLabels = 'cool warm hot'.split()
fields = 'm90 rho90Rvir rho90kpc'.split()
ylabels = [r'Mass Contained [M$_{\odot}$]',
            r'$r_{90}$ [Rvir]',
            r'$r_{90}$ [pkpc]']
            

expns = df.index
expns = [float(a.split('a')[-1])/1000. for a in expns]
redshifts = [1./a - 1. for a in expns]

colors = 'blue red green'.split()
markers = 'o ^ s'.split()

for tempLabel in tempLabels:

    fig1,ax1 = plt.subplots(1,1,figsize=(5,5))
    fig2,ax2 = plt.subplots(1,1,figsize=(5,5))
    fig3,ax3 = plt.subplots(1,1,figsize=(5,5))
    axes = [ax1,ax2,ax3]
    figs = [fig1,fig2,fig3]
    
    for denseLabel,color,mark in zip(denseLabels,colors,markers):
        for ax,field in zip(axes,fields):
            x = redshifts[:-5]
            y = df[tempLabel,denseLabel,field][:-5]
            ax.plot(x,y,color=color,marker=mark,label=denseLabel)

    for ax,ylab in zip(axes,ylabels):
        ax.set_xlabel('Redshift')
        ax.set_ylabel(ylab)
        ax.legend(loc='best')
        ax.invert_xaxis()
        if 'Mass' in ylab:
            ax.set_yscale('log')
            ax.set_ylim([10**6,10**10])
        elif 'kpc' in ylab:
            ax.set_ylim([0,250])
        else:
            ax.set_ylim([0,3])
        mkLines(ax)

    for fig,field in zip(figs,fields):
        s = 'filament_{0:s}_{1:s}.png'.format(field,tempLabel)
        fig.savefig(s,bbox_inches='tight',dpi=300)
        plt.close(fig)


for denseLabel in denseLabels:

    fig1,ax1 = plt.subplots(1,1,figsize=(5,5))
    fig2,ax2 = plt.subplots(1,1,figsize=(5,5))
    fig3,ax3 = plt.subplots(1,1,figsize=(5,5))
    axes = [ax1,ax2,ax3]
    figs = [fig1,fig2,fig3]
    
    for tempLabel,color,mark in zip(tempLabels,colors,markers):
        for ax,field in zip(axes,fields):
            x = redshifts[:-5]
            y = df[tempLabel,denseLabel,field][:-5]
            ax.plot(x,y,color=color,marker=mark,label=tempLabel)

    for ax,ylab in zip(axes,ylabels):
        ax.set_xlabel('Redshift')
        ax.set_ylabel(ylab)
        ax.legend(loc='best')
        ax.invert_xaxis()
        if 'Mass' in ylab:
            ax.set_yscale('log')
            ax.set_ylim([10**6,10**10])
        elif 'kpc' in ylab:
            ax.set_ylim([0,250])
        else:
            ax.set_ylim([0,3])
        mkLines(ax)

    for fig,field in zip(figs,fields):
        s = 'filament_{0:s}_{1:s}.png'.format(field,denseLabel)
        fig.savefig(s,bbox_inches='tight',dpi=300)
        plt.close(fig)


# Plot masses summed across all densities    
fig,ax = plt.subplots(1,1,figsize=(5,5))
for tempLabel,color,mark in zip(tempLabels,colors,markers):
    field = 'm90'
    x = redshifts[:-5]
    y = df.xs(field,level=2,axis=1)[tempLabel].sum(axis=1)[:-5]
    y = np.log10(y)
    ax.plot(x,y,color=color,marker=mark,label=tempLabel)
ax.set_xlabel('Redshift')
ax.set_ylabel(r'Log Mass [M$_{\odot}$]')
ax.legend(loc='best')
ax.invert_xaxis()
s = 'filament_totalMass_temperature.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)
     
# Plot masses summed across all temperatures
fig,ax = plt.subplots(1,1,figsize=(5,5))
for denseLabel,color,mark in zip(denseLabels,colors,markers):
    field = 'm90'
    x = redshifts[:-5]
    y = df.xs(field,level=2,axis=1).xs(denseLabel,
                level=1,axis=1).sum(axis=1)[:-5]
    y = np.log10(y)
    ax.plot(x,y,color=color,marker=mark,label=denseLabel)
ax.set_xlabel('Redshift')
ax.set_ylabel(r'Log Mass [M$_{\odot}$]')
ax.legend(loc='best')
ax.invert_xaxis()
s = 'filament_totalMass_density.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
    













