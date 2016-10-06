
from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import matplotlib.ticker as ticker

def fmt(x,pos):
    return '{0:.2f}'.format(x)

def mkHist(ax,x,y,z,cbarLabel,stat,binrange):

    numbins = 50

    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,z,
                                statistic=stat, bins=numbins,
                                range=binrange)

    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    if cbarLabel=='SNII':
        mesh = ax.pcolormesh(xedges,yedges,h,vmin=-7,vmax=-1)
    else:
        mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel('Rho [Rvir]')
    ax.set_ylabel('Z [Rvir]')
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True,
                        format=ticker.FuncFormatter(fmt))
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('{0:s} {1:s}'.format(stat,cbarLabel),rotation=270,fontsize=12)
    






pd.options.mode.chained_assignment = None  # default='warn'

baseloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

loTemps = [3.5,4.5,5.5]
hiTemps = [4.5,5.5,6.5]
tempLabels = ['cool','warm','hot']

limits = pd.read_csv('cloudLimits.csv')

expns0 = range(20,50)
expns = [i/100. for i in expns0]

minSNII = 2
maxSNII = -10

for i in range(len(loTemps)):
    print(tempLabels[i])
    loT = 10**loTemps[i]
    hiT = 10**hiTemps[i]
    for a in expns:
        
        print(a)

        with open(baseloc+rotname.format(a)) as f:
            f.readline()
            rvir = float(f.readline().split()[3])

        fname = baseloc+filename.format(a)
        df = pd.read_hdf(fname, 'data')


        nRange = limits[limits.expn==a]
        loN = 10**(nRange['loN'].values[0])
        hiN = 10**(nRange['hiN'].values[0])

        tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
        densInd = (df['density']>loN) & (df['density']<hiN)
        spacInd = (df['zRot']>0) & (df['theta']<80)

        cl = df[tempInd & densInd & spacInd]

        cl['xRotScale'] = cl['xRot']/rvir
        cl['yRotScale'] = cl['yRot']/rvir
        cl['zRotScale'] = cl['zRot']/rvir
        
        cl['rho'] = np.sqrt( cl['xRotScale']**2 + cl['yRotScale']**2 )

        plotting = 1
        if plotting == 1:
            fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(15,10))
            #ax1.scatter(cl['xRotScale'],cl['yRotScale'],marker='o',color='green',alpha=0.01)
            #ax2.scatter(cl['xRotScale'],cl['zRotScale'],marker='o',color='green',alpha=0.01)
            #ax3.scatter(cl['yRotScale'],cl['zRotScale'],marker='o',color='green',alpha=0.01)

            binrange = [[-3,3],[-3,3]]
            mkHist(ax1,cl['xRotScale'],cl['yRotScale'],cl['x'],'Count','count',binrange)    
            binrange = [[-3,3],[0,3]]
            mkHist(ax2,cl['xRotScale'],cl['zRotScale'],cl['x'],'Count','count',binrange)    
            binrange = [[-3,3],[0,3]]
            mkHist(ax3,cl['yRotScale'],cl['zRotScale'],cl['x'],'Count','count',binrange)    

            ax1.set_xlabel('x [Rvir]')
            ax2.set_xlabel('x [Rvir]')
            ax3.set_xlabel('y [Rvir]')

            ax1.set_ylabel('y [Rvir]')
            ax2.set_ylabel('z [Rvir]')
            ax3.set_xlabel('z [Rvir]')

            ax1.set_xlim([-3,3])
            ax2.set_xlim([-3,3])
            ax3.set_xlim([-3,3])

            ax1.set_ylim([-3,3])
            ax2.set_ylim([0,3])
            ax3.set_ylim([0,3])

            ax2.set_title('a = {0:.3f}, z = {1:.3f}'.format(a,1./a-1))
            
            # Plot Historgrams
            binrange = [[0,3],[0,3]]
            mkHist(ax4,cl['rho'],cl['zRotScale'],cl['SNII'],'SNII','mean',binrange)
            mkHist(ax5,cl['rho'],cl['zRotScale'],cl['density'],'nH','mean',binrange)
            mkHist(ax6,cl['rho'],cl['zRotScale'],cl['temperature'],'Temp','mean',binrange)

            fig.tight_layout()

            s = 'rotation_test_a{0:.3f}_{1:s}Gas.png'.format(a,tempLabels[i])
            fig.savefig(s,bbox_inches='tight',dpi=300)
            plt.close()
            
        currentMinSNII = np.log10(cl['SNII'].min())
        currentMaxSNII = np.log10(cl['SNII'].max())

        if currentMinSNII < minSNII:
            minSNII = currentMinSNII
        if currentMaxSNII > maxSNII:
            maxSNII = currentMaxSNII


    print('SNII Range = {0:.3f} - {1:.3f}'.format(minSNII, maxSNII))
    


    



