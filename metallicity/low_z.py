
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

def spread(x):
    return max(x) - min(x)


baseloc = '/mnt/cluster/abs/cgm/vela2b/'

galnum = range(21,30)
expn = ['0.490']*len(galnum)
expn[galnum.index(24)] = '0.450'

fig1, ((ax11,ax21,ax31),(ax41,ax51,ax61),(ax71,ax81,ax91)) = plt.subplots(3,3,figsize=(15,15))
axes1 = (ax11, ax21, ax31, ax41, ax51, ax61, ax71, ax81, ax91)
fig2, ((ax12,ax22,ax32),(ax42,ax52,ax62),(ax72,ax82,ax92)) = plt.subplots(3,3,figsize=(15,15))
axes2 = (ax12, ax22, ax32, ax42, ax52, ax62, ax72, ax82, ax92)

lowT, highT = 4, 4.5
lown, highn = -5, -4
coldStats = '{0:s}vela2b_lowZgas_stats.out'.format(baseloc)
f = open(coldStats, 'w')
header = '# Stats of Metalicty of cells within phase cut\n # Cut = {0:.1f} < log(T) < {1:.1f}\t{2:.1f} < log(nH) < {3:.1f}\n'.format(lowT,highT,lown,highn) 
f.write(header)
header = 'Galnum\tnObs\tMin\tMax\tMean\tVar\tSkew\tKurt\n'
f.write(header)
s = '{0:d}\t{1:d}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{5:.6f}\t{6:.6f}\t{7:.6f}\n'

for galID, a, ax1, ax2 in zip(galnum, expn, axes1, axes2):

    boxfile = '{0:s}vela{1:d}/a{2:s}/vela2b-{1:d}_GZa{2:s}.tcdat.h5'.format(baseloc,galID,a)
    print(boxfile)

    df = pd.read_hdf(boxfile, 'data')
    n = np.log10(df['nH'])
    t = np.log10(df['T'])
    z = np.log10(df['Zcell'])
    ids = np.arange(1,len(n))
    
    index = (t>4) & (t<4.5) & (n>-5) & (n<-4)
    
    cold_z = z[index]
    hot_z = z[np.invert(index)]
    c = stats.describe(cold_z)
    f.write(s.format(galID,c[0],c[1][0],c[1][1],c[2],c[3],c[4],c[5]))

    ax1.hist(cold_z, bins=20, log=True, histtype='step', label='Cloud')
    ax1.hist(hot_z, bins=20, log=True, histtype='step', label='Rest')
    ax1.legend(frameon=False, loc='upper left')
    ax1.set_xlabel('Metallicity')
    ax1.set_title('Vela2b-{0:d}'.format(galID))

    binrange = [[np.floor(n.min()),np.ceil(n.max())],[np.floor(t.min()),np.ceil(t.max())]]
    binrange = [[-9, 3],[2,8]]
    h, xedges, yedges, binnumber = stats.binned_statistic_2d(n,t,z, statistic=spread, bins=50, range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0, h)
     
    mesh = ax2.pcolormesh(xedges, yedges, h, cmap='viridis')
    ax2.set_xlim(binrange[0])
    ax2.set_ylim(binrange[1])
    print( '\t',binrange[1] )
    ax2.set_xlabel('Density')
    ax2.set_ylabel('Temperature')
    ax2.set_title('Vela2b-{0:d}'.format(galID))
    cbar = fig2.colorbar(mesh, ax=ax2, use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Z', rotation=270, fontsize=12)

fig1.tight_layout()
fig1.savefig('vela2b_metallicityHist.png', bbox_inches='tight', dpi=300)
    
fig2.tight_layout()
fig2.savefig('vela2b_phase_Zspread.png', bbox_inches='tight', dpi=300)
    

