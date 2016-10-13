
'''
Calculates SFR for vela2b-27
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def getMass(a):

    dataloc = '/mnt/cluster/abs/Simulations/vela2.1/VELA27/output/ana/'
    dataloc = '/home/jacob/research/velas/vela2b/vela27/subhalos/'
    filename = 'halos_{0:.3f}.txt'
    fname = dataloc + filename.format(a)

    with open(fname) as f:
        f.readline()
        f.readline()
        l = f.readline().split()
        mstar = float(l[13])
        mvir = float(l[7])
        rvir = float(l[8])
        mgas = float(l[12])

    return mstar,mvir,rvir,mgas


times = pd.read_csv('lookback_time.csv',index_col='a')

# Initial mass
a = 0.200
mstar1,mvir1,rvir1,mgas1 = getMass(a)
timestep = times.ix[a]['Time until next [Myr]']*1e6

expns = [i/100. for i in range(21,50)]

header = ['a','mstar','mvir','sfr','ssfr']
data = np.zeros(len(header))

for a in expns:

    mstar2,mvir2,rvir2,mgas2 = getMass(a)

    print(mstar1,mstar2)
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
df.to_csv('vela2b-27_sfr.csv',index=False)



fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.plot(df['a'],df['sfr'])
ax2.plot(df['a'],df['ssfr'])
ax1.set_ylabel('SFR')
ax2.set_ylabel('sSFR')
ax1.set_xlabel('a')
ax2.set_xlabel('a')
# Merger event
ax1.vlines(0.29,df['sfr'].min(),df['sfr'].max(),colors='red',linestyle='dashed',linewidth=2)
ax1.vlines(0.32,df['sfr'].min(),df['sfr'].max(),colors='red',linestyle='dashed',linewidth=2)
ax2.vlines(0.29,df['ssfr'].min(),df['ssfr'].max(),colors='red',linestyle='dashed',linewidth=2)
ax2.vlines(0.32,df['ssfr'].min(),df['ssfr'].max(),colors='red',linestyle='dashed',linewidth=2)
# Dispersion event
ax1.vlines(0.35,df['sfr'].min(),df['sfr'].max(),colors='green',linestyle='dashed',linewidth=2)
ax1.vlines(0.45,df['sfr'].min(),df['sfr'].max(),colors='green',linestyle='dashed',linewidth=2)
ax2.vlines(0.35,df['ssfr'].min(),df['ssfr'].max(),colors='green',linestyle='dashed',linewidth=2)
ax2.vlines(0.45,df['ssfr'].min(),df['ssfr'].max(),colors='green',linestyle='dashed',linewidth=2)
fig.tight_layout()
fig.savefig('vela2b-27_sf.png',bbox_inches='tight',dpi=300)







