
'''
Plots the evolution of the phase diagram of vela2b
'''


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


baseloc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:.3f}/'
basefile = '{0:s}vela2b-{1:d}_GZa{2:.3f}.h5'

numbins = 50
binrange = [[-10,2],[2,8]]

galNums = range(21,30)
expns = np.arange(0.200, 0.500, 0.01)

for galNum in galNums:

    print('Galnum = {0:d}'.format(galNum))

    for a in expns:
        
        print('\ta = {0:.3f}'.format(a))
        redshift = 1./a - 1

        loc = baseloc.format(galNum,a)
        fname = basefile.format(loc,galNum,a)

        try:
            df = pd.read_hdf(fname, 'data')
        except IOError:
            continue

        
        # Bin
        H, xedges, yedges = np.histogram2d(np.log10(df['density']), np.log10(df['temperature']), 
                                            bins=numbins, range=binrange)
        H = np.rot90(H)
        H = np.flipud(H)
        H = np.ma.masked_where(H==0,H)
        H = np.log10(H)

        fig, ax = plt.subplots(1,1,figsize=(5,5))
        mesh = ax.pcolormesh(xedges,yedges,H)
        ax.set_xlim(binrange[0])
        ax.set_ylim(binrange[1])
        ax.set_xlabel('Density')
        ax.set_ylabel('Temperature')
        ax.set_title('z = {0:.2f}'.format(redshift))
        
        cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel('Counts', rotation=270, fontsize=12)


        s = 'vela2b-{0:d}_a{1:.3f}_phase.png'.format(galNum,a)
        fig.savefig(s,bbox_inches='tight',dpi=300)




