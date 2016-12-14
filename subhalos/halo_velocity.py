
'''
Plots the velocity of the host halo and the largest subhalo over time
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


dataloc = '/home/jacob/research/velas/vela2b/vela27/subhalos/'
filename = 'halos_{0:.3f}.txt'

expns = np.arange(0.200,0.550,0.01)


hs, ss, ds, vr, d = [], [], [], [], []

for a in expns:

    fname = dataloc+filename.format(a)

    data = np.loadtxt(fname,skiprows=2)

    vx = data[:,4]
    vy = data[:,5]
    vz = data[:,6]
    mvir = data[:,7]
    rvir = data[:,8]
    xc = data[:,15]
    yc = data[:,16]
    zc = data[:,17]
    
    vxDiff = vx[1]-vx[0]
    vyDiff = vy[1]-vy[0]
    vzDiff = vz[1]-vz[0]

    hostSpeed = np.sqrt(vx[0]**2 + vy[0]**2 + vz[0]**2)
    subSpeed = np.sqrt(vx[1]**2 + vy[1]**2 + vz[1]**2)
    diffSpeed = np.sqrt(vxDiff**2 + vyDiff**2 + vzDiff**2)

    distance = np.sqrt(xc[1]**2 + yc[1]**2 + zc[1]**2)
    vrad = (vxDiff*xc[1] + vyDiff*yc[1] + vzDiff*zc[1]) / distance
    
    hs.append(hostSpeed)
    ss.append(subSpeed)
    ds.append(diffSpeed)
    vr.append(vrad)
    d.append(distance)

    
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))

ax1.plot(expns,hs,label='host')
ax1.plot(expns,ss,label='subhalo')
ax1.plot(expns,ds,label='relative')
ax1.plot(expns,vr,label='radial')
ax2.plot(expns,d,label='distance')

for ax in [ax1,ax2]:
    ax.set_xlabel('Expansion Parameter')
    ymin,ymax = ax.get_ylim()
    ax.vlines(0.29,ymin,ymax,linestyle='dashed')
    ax.vlines(0.32,ymin,ymax,linestyle='dashed')
    ax.vlines(0.35,ymin,ymax,linestyle='dashed')

ax1.set_ylabel('Velocity [km/s]')
ax2.set_ylabel('Distance [kpc]')

ax1.legend(loc='best')

fig.tight_layout()
fig.savefig('haloVelocities.png',bbox_inches='tight',dpi=300)





