
import pandas as pd
import matplotlib.pyplot as plt

eps = 50
minPart = 1000

galNums = range(21,30)

fname = 'vela2b_onion_DBSCAN_eps{0:d}_min{1:d}_velocity.h5'.format(eps,minPart)
df = pd.read_hdf(fname, 'data')

for galNum in galNums:
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))

    d = df[(df['GalNum']==galNum).values]
    ax1.plot(d['Cluster'], d['vr_min'], 'kx', label='min')
    ax1.plot(d['Cluster'], d['vr_max'], 'k+', label='max')
    ax1.plot(d['Cluster'], d['vr_mean'], 'ko', label='mean')
    ax1.plot(d['Cluster'], d['vr_std']+d['vr_mean'], 'k_', label='std')
    ax1.plot(d['Cluster'], d['vr_std']-d['vr_mean'], 'k_')

    ax2.plot(d['Cluster'], d['r_min'], 'kx', label='min')
    ax2.plot(d['Cluster'], d['r_max'], 'k+', label='max')
    ax2.plot(d['Cluster'], d['r_mean'], 'ko', label='mean')
    ax2.plot(d['Cluster'], d['r_std']+d['vr_mean'], 'k_', label='std')
    ax2.plot(d['Cluster'], d['r_std']-d['vr_mean'], 'k_')

    ax1.set_xlabel('Cluster Number')
    ax1.set_ylabel('Radial Velocity [km/s]')
    ax1.legend(frameon=False, fontsize='small')
    ax2.set_xlabel('Cluster Number')
    ax2.set_ylabel('Radial Location [kpc]')
    ax2.legend(frameon=False)

    fig.tight_layout()
    s = 'vela2b-{2:d}_onion_DBSCAN_eps{0:d}_min{1:d}_velocity.png'.format(eps,minPart,galNum)
    fig.savefig(s, bbox_inches='tight')
    






