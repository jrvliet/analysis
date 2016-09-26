
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import gc
import statsmodels as sm
import sys


mH = 1.6737e-24
pc2cm = 3.086e18
mSun = 1.989e33

vmin, vmax = -8, -1

galNums = range(21,30)
galNums = range(27,30)
galNums = [27]
expns = np.arange(0.200,0.500,0.01) 
#expns = ['{0:.3f}'.format(i) for i in a]

for galNum in galNums:
    for a in expns:
        print a

        dataloc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:.3f}/'.format(galNum,a)
        dataloc = '/home/jacob/research/velas/vela2b/vela{0:d}/a{1:.3f}/'.format(galNum,a)
        boxfile = '{0:s}vela2b-{1:d}_GZa{2:.3f}.h5'.format(dataloc,galNum,a)
        rotmatFile = '{0:s}rotmat_a{1:.3f}.txt'.format(dataloc,a)
        try:
            f = open(rotmatFile, 'r')
        except:
            continue

        f.readline()
        l = f.readline().split()
        rvir = float(l[3])
        f.close()
            
        d = pd.read_hdf(boxfile, 'data')

        dense = np.arange(-5,-2,0.25)
        loT, hiT = 10**3, 10**4.5
        loN = 10**-6

        baseInds = ( (d['temperature']<=hiT) & (d['temperature']>=loT) & 
                     (d['x']<0) & (d['z']>0) & (np.abs(d['y'])<300))
        #baseInds = ( (d['temperature']<=hiT) & (d['temperature']>=loT) )

        fig = plt.figure(figsize=(20,15))
        axes = []
        for i,n in enumerate(dense):
            print i
            if i>=1:
                loN = 10**(dense[i-1])
            hiN = 10**n
            cloudInds = (baseInds & (d['density']<=hiN) & (d['density']>=loN))
            print 'log(loN) = {0:.2f}\tlog(hiN) = {1:.2f}\tsum(cloudInds) = {2:d}'.format(np.log10(loN),
                    np.log10(hiN),sum(cloudInds))
            cloud = d[cloudInds]

            # Calculate mass weighted metallicity mean
            cloud['mass'] = cloud['density']*mH* (cloud['cell_size']*pc2cm)**3 / mSun

            cloud['weight'] = cloud['mass']*cloud['SNII']
            cloud['meanZ'] = np.log10(cloud['weight'].sum() / cloud['mass'].sum())
            print cloud['meanZ'].min(), cloud['meanZ'].max()


            ax = fig.add_subplot(3, 4, i+1, projection='3d')
            axes.append(ax)
            im = ax.scatter(cloud['x']/rvir, cloud['y']/rvir, cloud['z']/rvir, c=cloud['meanZ'], 
                        vmin=vmin, vmax=vmax, marker='o', alpha=0.01)
            ax.view_init(elev=90, azim=0)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            ax.set_xlim([-3,0])
            ax.set_ylim([-3,3])
            ax.set_zlim([0,3])
            ax.view_init(elev=30, azim=0)
            ax.set_title('{0:.2f} $<$ log(nH) $<$ {1:.2f}'.format(np.log10(loN),np.log10(hiN)))



        z = 1.0/float(a) - 1
        fig.suptitle('a = {0:.2f}, Rvir = {1:.1f}'.format(a,rvir))
        fig.tight_layout()
    
        norm = mp.colors.Normalize(vmin=vmin,vmax=vmax)
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85,0.05,0.02,0.9])
        mp.colorbar.ColorbarBase(cbar_ax, cmap='viridis', norm=norm)
        #fig.colorbar(im,cax=cbar_ax)
        #cax, kw = mp.colorbar.make_axes(axes)
        #cbar = fig.colorbar(im,cax=cax, **kw)
        #cbar.set_clim(-8,-1)

        s = './density_plots/vela2b-{0:d}_a{1:.3f}_cloudSize.png'.format(galNum,a)
        fig.savefig( s, bbox_inches='tight', pdi=300)
        



