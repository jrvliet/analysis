


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plotStat(df,stat,tempLabel):
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    for name in df.columns:
        if stat in name:
            ax.plot(df['a'],np.log10(df[name]),label=name.strip(stat)+'Z')
            print('{0:s}\t{1:s}\t{2:.3f}-{3:.3f}'.format(name,tempLabel,
                    np.log10(df[name].min()),np.log10(df[name].max())))

    ax.set_xlabel('Expansion Parameter')
    
    ylab = r'log($P_{{{0:s}}})$'.format(stat.lower())
    if 'mean' in stat.lower():
        ylab = r'log($\bar{P})$'
        ax.set_ylim([-16,-13])
    elif 'std' in stat.lower():
        ylab = r'log($\sigma_{P}$)'
        ax.set_ylim([-17,-12])
    ax.set_ylabel(ylab)
    #ax.set_yscale('log')
    ax.legend(loc='best')
    ymin,ymax = ax.get_ylim()
    ax.vlines(0.29,ymin,ymax,linestyle='dashed')
    ax.vlines(0.32,ymin,ymax,linestyle='dashed')
    ax.vlines(0.35,ymin,ymax,linestyle='dashed')
    s = 'pressureStats_{0:s}_{1:s}.png'.format(stat.lower(),tempLabel)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)




dataloc = '/home/jacob/research/code/analysis/inflows/'
filename = 'pressureStats_noCGM_{0:s}.h5'

tempLabels = 'cool warm hot'.split()
for tempLabel in tempLabels:
    fname = dataloc+filename.format(tempLabel)

    df = pd.read_hdf(fname,'data')
    df['a'] = np.arange(0.200,0.550,0.01)

    plotStat(df,'Mean',tempLabel)
    plotStat(df,'Std',tempLabel)


