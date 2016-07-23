
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


baseloc = '/mnt/cluster/abs/cgm/vela2b/'
galnum = range(21,30)
expn = ['0.490']*len(galnum)
expn[galnum.index(24)] = '0.450'

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(15,15))
axes = (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9)

for galID, a, ax in zip(galnum, expn, axes):

    loc = '{0:s}vela{1:d}/a{2:s}/'.format(baseloc,galID,a)
    fname = '{0:s}vela2b-{1:d}_GZa{2:s}.h5'.format(loc,galID,a)

    dfull = pd.read_hdf(fname, 'data')
    zCut = dfull.quantile(0.1)['SNII']
    d = dfull[(dfull.SNII<zCut)]


    r = np.sqrt(d['x']**2 + d['y']**2 + d['z']**2)
    rfull = np.sqrt(dfull['x']**2 + dfull['y']**2 + dfull['z']**2)
        
    vrad = (d['vx']*d['x'] + d['vy']*d['y'] + d['vz']*d['z'] )/r
    vradfull = (dfull['vx']*dfull['x'] + dfull['vy']*dfull['y'] + dfull['vz']*dfull['z'] )/rfull

    ax.hist(vrad, bins=20, histtype='step', log=True, color='blue', label='Z$<${0:.1f}'.format(np.log10(zCut)))
    ax.hist(vradfull, bins=20, histtype='step', log=True, color='green', label='Full Box')
    ax.legend(frameon=False)
    ax.set_xlabel('Radial Velocity [km/s]')
    ax.set_title('Vela2b-{0:d}'.format(galID))

fig.tight_layout()
fig.savefig('vela2b_vrad_lowZcut.png', bbox_inches='tight', dpi=300)

    


