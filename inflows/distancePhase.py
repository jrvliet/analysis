
'''
Plots the phase diagram of the filament region. 
Color each cell is the mean galactocentric distance
'''


from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
pd.options.mode.chained_assignment = None
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import sys


def fmt(x,pos):
    return '{0:.1f}'.format(x)


def mkHist(ax,x,y,z,stat,xlabel,ylabel,lox,hix,loy,hiy,field):
    numbins = 50
    #stat = 'count'
    binrange = [[lox,hix],[loy,hiy]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,z,
                              statistic=stat, bins=numbins,
                                 range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0

    #vmax = -1
    #vmin = -6
    #if vmax<h.max():
    #    print('Exceeds\t',h.max())
    h = np.ma.masked_where(h==0,h)
#    h = np.log10(h)
    print(h.max(),np.percentile(h,10))

    mesh = ax.pcolormesh(xedges,yedges,h)
    #mesh = ax.pcolormesh(xedges,yedges,h,vmin=0,vmax=vmax)
    #mesh = ax.pcolormesh(xedges,yedges,h,vmin=vmin,vmax=vmax,cmap='jet')
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True,
                       format=ticker.FuncFormatter(fmt) )
    cbar.ax.get_yaxis().labelpad = 20
    if 'Mod' in field:
        clabel = 'Mean Distance [Rvir]'
    else:
        clabel = 'Mean Distance [pkpc]'
    cbar.ax.set_ylabel(clabel,rotation=270,fontsize=12)
    return h,xedges,yedges,binnumber



loc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
loc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'
expns = np.arange(0.210,0.540,0.01)


sfrLoc = '/home/jacob/research/code/analysis/sfr/'
sfrFile = 'vela2b-27_sfr.csv'
sfr = pd.read_csv(sfrLoc+sfrFile)
sfr['redshift'] = 1./sfr['a'] - 1

# In[4]:

loT,hiT = 3.5,6.5
loN,hiN = -5.5,-2.5

header = 'a nMax tMax'.split()
maxMetal = np.zeros((len(expns),len(header)))
maxMetal = pd.DataFrame(maxMetal,columns=header)

# In[18]:

for i,a in enumerate(expns):
    print(a)
    redshift = 1./a - 1
    
    # Get rvir
    rname = loc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    # Read in gas box
    fname = loc+filename.format(a)
    df = pd.read_hdf(fname,'data')

    # Calclate distance
    df['r'] = np.sqrt(df['xRot']**2 + df['yRot']**2 + df['zRot']**2)
    df['rMod'] = df['r']/rvir
    
    # Select out filament
    tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    densInds = (df['density']>10**loN) & (df['density']<10**hiN)
    spacInds = (df['theta']<80)
    fil = df[tempInds & densInds & spacInds]

    fields = 'r rMod'.split()
    units = 'pkpc rvir'.split()
    for field,unit in zip(fields,units):
        fig = plt.figure(figsize=(5,7))
        gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
        gs.update(hspace=0.25)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])

        # Phase Plot
        h,xedges,yedges,binnumber = mkHist(ax1,np.log10(fil['density']),
                                    np.log10(fil['temperature']),fil[field],
                                    'mean','Density','Temperature',loN,hiN,loT,hiT,field)
        ax1.set_title('{0:.3f}'.format(redshift))

        # SFR Plot
        ax2.plot(sfr['redshift'],sfr['sfr'],color='black')
        ax2.scatter(sfr['redshift'].ix[i],sfr['sfr'].ix[i],color='red',marker='s')
        ax2.invert_xaxis()
        ax2.set_xlabel('Redshift')
        ax2.set_ylabel(r'SFR [M$_{\odor}$ yr$^{-1}$]')
        ymin,ymax = ax2.get_ylim()
        
        # Save and close
        s = 'distancePhase_a{0:.3f}_noCGM_{1:s}.png'.format(a,unit)
        fig.savefig(s,bbox_inches='tight',dpi=300)
        plt.close(fig)



