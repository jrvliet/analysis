
from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
from numpy import log10
import sys


baseloc = '/mnt/cluster/abs/cgm/vela2b/'
galnum = range(21,30)
expn = ['0.490']*len(galnum)
expn[galnum.index(24)] = '0.450'

fout = 'vela2b_metallicity_stats.dat'
f = open(fout, 'w')
header = '\t\tSNII\t\t\t\t\tZ\nGalID\t\tMin\tMax\tMean\tStd\t10%\tMin\tMax\tMean\tStd\t10%\n'
f.write(header)
s = '{0:s}\t{1:.3f}\t{2:.3f}\t{3:.3f}\t{4:.3f}\t{5:.3f}\t{6:.3f}\t{7:.3f}\t{8:.3f}\t{9:.3f}\t{10:.3f}\n'


fig, ((ax1, ax2, ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(15,15))
axes = (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9)

for galID,a,ax in zip(galnum,expn,axes):

    loc = '{0:s}vela{1:d}/a{2:s}/'.format(baseloc,galID,a)

    boxfile = '{0:s}vela2b-{1:d}_GZa{2:s}.h5'.format(loc,galID,a)
    tcdatfile = '{0:s}vela2b-{1:d}_GZa{2:s}.tcdat.h5'.format(loc,galID,a)

    print(boxfile)

    box = pd.read_hdf(boxfile, 'data')
    cool = pd.read_hdf(tcdatfile, 'data')

    snIIStats = log10(box['SNII'].describe())
    ZStats = log10(cool['Zcell'].describe())
    snIIpercent = log10(box['SNII'].quantile(0.1))
    Zpercent = log10(cool['Zcell'].quantile(0.1))
    galname = 'vela2b-{0:d}'.format(galID)
    f.write(s.format(galname,snIIStats[3],snIIStats[7],snIIStats[1],snIIStats[2],snIIpercent,
                     ZStats[3],ZStats[7],ZStats[1],ZStats[2],Zpercent))

    ax.loglog(box['SNII'], cool['Zcell'], 'k.')
    ax.set_xlabel('SNII MF')
    ax.set_ylabel('Zcell')
    ax.set_title('Vela2b-{0:d}'.format(galID))

fig.tight_layout()
fig.savefig('metallicity_comparison.png', bbox_inches='tight', dpi=300)

f.close()

