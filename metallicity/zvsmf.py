
from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt


baseloc = '/mnt/cluster/abs/cgm/vela2b/'
galnum = range(21,30)
expn = ['0.490']*len(galnum)
expn[galnum.index(24)] = '0.450'


fig, ((ax1, ax2, ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(15,15))
axes = (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9)

for galID,a,ax in zip(galnum,expn,axes):

    loc = '{0:s}vela{1:d}/a{2:s}/'.format(baseloc,galID,a)

    boxfile = '{0:s}vela2b-{1:d}_GZa{2:s}.h5'.format(loc,galID,a)
    tcdatfile = '{0:s}vela2b-{1:d}_GZa{2:s}.tcdat.h5'.format(loc,galID,a)

    print(boxfile)

    box = pd.read_hdf(boxfile, 'data')
    cool = pd.read_hdf(tcdatfile, 'data')

    ax.loglog(box['SNII'], cool['Zcell'], 'k.')
    ax.set_xlabel('SNII MF')
    ax.set_ylabel('Zcell')
    ax.set_title('Vela2b-{0:d}'.format(galID))

fig.tight_layout()
fig.savefig('metallicity_comparison.png', bbox_inches='tight', dpi=300)


