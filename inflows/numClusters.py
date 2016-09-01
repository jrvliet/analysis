
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

dataloc = '/home/jacob/research/code/analysis/inflows/stats/'
basename = '{0:s}/vela2b-{1:d}_onion_DBSCAN_stats.h5'

galNums = range(21,30)

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(20,15))
axes = (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9)

for galNum, ax in zip(galNums,axes):

    print('Galnum = {0:d}'.format(galNum))
    df = pd.read_hdf(basename.format(dataloc,galNum), 'data')
    
    nEPS = len(df['eps'].unique())-1
    nMin = len(df['minPart'].unique())-1
    n = np.zeros((nEPS,nMin))
    eps = df['eps'].unique()[1:]
    minP = df['minPart'].unique()[1:]
    for i in range(len(eps)):
        for j in range(len(minP)):
            index = ((df['eps']==eps[i]) & (df['minPart']==minP[j]))
            n[i,j] = sum(index)

    im = ax.imshow(n, cmap='viridis',origin='lower',aspect='auto',
                    extent=(minP[0],minP[-1],eps[0],eps[-1]))
    ax.set_xlabel('Min Particles')
    ax.set_ylabel('EPS')
    ax.set_title('Vela {0:d}'.format(galNum))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.get_yaxis().labelpad = 12
    cbar.ax.set_ylabel('Number of clusters',rotation=270)
    #divider = make_axes_locatable(vars()[ax])
    #vars()["c"] = divider.append_axes("right", size = "5%", pad = 0.05)

#fig.tight_layout()
s = '{0:s}/vela2b_dbscan_stats.png'.format(dataloc)
fig.savefig(s, bbox_inches='tight', dpi=300)



