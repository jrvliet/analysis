
'''
Calculates SFR for vela2b-27
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys


def getMass(a,galNum):

    dataloc = '/home/jacob/research/velas/vela2b/vela27/subhalos/'
    dataloc = '/mnt/cluster/abs/Simulations/vela2.1/VELA{0:d}/output/ana/'
    filename = 'halos_{0:.3f}.txt'
    fname = dataloc.format(galNum) + filename.format(a)

    try:
        with open(fname) as f:
            f.readline()
            f.readline()
            l = f.readline().split()
            mstar = float(l[13])
            mvir = float(l[7])
            rvir = float(l[8])
            mgas = float(l[12])

        return mstar,mvir,rvir,mgas
    except IOError:
        return 0,0,0,0


def mkLines(ax):
    ymin,ymax = ax.get_ylim()
    points = [2.448,2.125,1.857]
    for point in points:
        ax.vlines(point,ymin,ymax,linestyle='dashed',linewidth=1,color='black')
    ax.set_ylim([ymin,ymax]) 
    

times = pd.read_csv('lookback_time.csv',index_col='a')

galNums = range(21,30)
    
for galNum in galNums:

    print(galNum)
    # Initial mass
    a = 0.200
    mstar1,mvir1,rvir1,mgas1 = getMass(a,galNum)
    timestep = times.ix[a]['Time until next [Myr]']*1e6

    expns = [i/100. for i in range(21,55)]

    header = ['a','mstar','mvir','sfr','ssfr']
    data = np.zeros(len(header))

    for a in expns:

        mstar2,mvir2,rvir2,mgas2 = getMass(a,galNum)
        if mstar2==0:
            break

        deltaMs = mstar2-mstar1
        sfr = deltaMs/timestep
        ssfr = sfr/mstar2

        
        d = np.zeros(len(header))
        d[0] = a
        d[1] = np.log10(mstar2)
        d[2] = np.log10(mvir2)
        d[3] = sfr
        d[4] = np.log10(ssfr)
        
        timestep = times.ix[a]['Time until next [Myr]']*1e6
        mstar1,mvir1,rvir1,mgas1 = mstar2,mvir2,rvir2,mgas2 

        data = np.vstack((data,d))

    data = np.delete(data,(0),axis=0)
    df = pd.DataFrame(data,columns=header)
    outname = 'vela2b-{0:d}_sfr.csv'.format(galNum)
    df.to_csv(outname,index=False)

sys.exit()
df['redshift'] = 1./df['a'] - 1.

# Plot SFR
fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(df['redshift'],df['sfr'],color='black')
ax.set_ylabel('SFR [M$_{\odot}$ yr$^{-1}$]')
ax.set_xlabel('Redshift')
ax.invert_xaxis()
mkLines(ax)
s = 'vela2b-27_sfr.pdf'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)

# Plot sSFR
fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(df['redshift'],df['ssfr'],color='black')
ax.set_ylabel('sSFR [yr$^{-1}$]')
ax.set_xlabel('Redshift')
ax.invert_xaxis()
mkLines(ax)
s = 'vela2b-27_ssfr.pdf'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)


