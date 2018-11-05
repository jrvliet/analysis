
import pandas as pd
import numpy as np
import itertools as it
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.stats as st



def mkHist(ax,x,y,z,stat,xlabel,ylabel,cbarax,numbins=100,smooth=False):
    #stat = 'count'
    limit = 325
    limit = 150
    binrange = [[-limit,limit],[-limit,limit]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,z,
                              statistic=stat, bins=numbins,
                                 range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    #h[np.isnan(h)] = np.nanmin(h)
    h[np.isnan(h)] = 0.
    h = np.ma.masked_where(h==0,h)
    #h = np.log10(h)
 
    mesh = ax.pcolormesh(xedges,yedges,h,cmap='viridis')#,vmin=vmin,vmax=vmax)
    if smooth:
        f = interp2d(x, y, z, kind='cubic')

        data1 = f(xnew,ynew)
        Xn, Yn = np.meshgrid(xedges, ynew)
        plt.subplot(3, 2, 5)
        plt.pcolormesh(Xn, Yn, data1, cmap='RdBu')
        
        
    mesh = ax.pcolormesh(xedges,yedges,h,cmap='jet')
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlabel,fontsize='x-large')
    ax.set_ylabel(ylabel,fontsize='x-large')
    
    print(h.max(),h.min())
    #return h,xedges,yedges,binnumber
    return mesh


# ## New file

# In[ ]:


loc = '/lustre/projects/p089_swin/jvander/analysis/'
filename = 'vela2b-{0:d}.a0.490.i90_OVIcellsGalFrame_abscells.h5'
galNums = range(21,30)
for galNum in galNums:
    fname = loc+filename.format(galNum)
    try:
        df = pd.read_hdf(fname,'data')
    except IOError:
        continue
    print(galNum)
    xlim = max([df['xGal'].max(),np.abs(df['xGal'].min())])
    ylim = max([df['yGal'].max(),np.abs(df['yGal'].min())])
    zlim = max([df['zGal'].max(),np.abs(df['zGal'].min())])

    df['cell_cm'] = df['cell_size']*3.0857e21

    df['logN'] = np.log10(df['cell_cm']*df['nH'].apply(lambda x: pow(10.,x)))
    df['logNOVI'] = np.log10(df['cell_cm']*df['ion_density'].apply(lambda x: pow(10.,x)))

    df[['logN','logNOVI']].describe()

    directions = [list('xy'),list('xz'),list('yz')]
    directions = [['xGal','yGal'], ['xGal','zGal'], ['yGal','zGal']]
    stats = 'median mean max'.split()
    for stat in stats:
        print(stat)
        fig,axes = plt.subplots(1,3,figsize=(15,5))
        cbarax = fig.add_axes([.95,.12,.03,.76])
        for ax,direction in zip(axes,directions):
            x = direction[0]
            y = direction[1]
            mesh = mkHist(ax,df[x],df[y],df['logNOVI'],stat,
                          x,y,cbarax,numbins=100)

        cbar = plt.colorbar(mesh,cax=cbarax,use_gridspec=True)
                            #format=ticker.FuncFormatter(fmt))
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel('logNOVI',rotation=270,fontsize='x-large')
        fig.subplots_adjust(wspace=.25)
        s = ('logNOVI_vela{0:d}_projection'
             '_{1:s}_galFrame_absCells.png'.format(galNum,stat))
        fig.savefig(s,dpi=300,bbox_inches='tight')
        plt.close(fig)

    directions = [list('xy'),list('xz'),list('yz')]
    directions = [['xGal','yGal'], ['xGal','zGal'], ['yGal','zGal']]
    stats = 'median mean max'.split()
    fields = 'vxGal vyGal vzGal vr'.split()
    for stat in stats:
        print(stat)
        fig,axes = plt.subplots(4,3,figsize=(15,15))
        cbarax = fig.add_axes([.95,.12,.03,.76])
        for i,direction in enumerate(directions):
            for j,field in enumerate(fields):
                x = direction[0]
                y = direction[1]
                ax = axes[j,i]
                mesh = mkHist(ax,df[x],df[y],df[field],stat,
                              x,y,cbarax,numbins=100)
                ax.set_title('{0:s}'.format(field))
        cbar = plt.colorbar(mesh,cax=cbarax,use_gridspec=True)
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel('logNOVI',rotation=270,fontsize='x-large')
        fig.subplots_adjust(wspace=.25,hspace=.35)
        s = 'velocity_vela{0:d}_projection_{1:s}_galFrame_absCells.png'.format(galNum,stat)
        fig.savefig(s,dpi=300,bbox_inches='tight')

        plt.close(fig)

    directions = [list('xy'),list('xz'),list('yz')]
    directions = [['xGal','yGal'], ['xGal','zGal'], ['yGal','zGal']]
    stats = 'median mean max'.split()
    fields = 'vr vtheta vphi'.split()
    rvir = df['rvir'].mean()
    d = df[df['r']<rvir]
    for stat in stats:
        print(stat)
        fig,axes = plt.subplots(3,3,figsize=(15,15))
        cbarax = fig.add_axes([.95,.12,.03,.76])
        for i,direction in enumerate(directions):
            for j,field in enumerate(fields):
                x = direction[0]
                y = direction[1]
                ax = axes[j,i]
                mesh = mkHist(ax,df[x],df[y],df[field],stat,
                              x,y,cbarax,numbins=100)
                ax.set_title('{0:s}'.format(field))
        cbar = plt.colorbar(mesh,cax=cbarax,use_gridspec=True)
                            #format=ticker.FuncFormatter(fmt))
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel('logNOVI',rotation=270,fontsize='x-large')
        fig.subplots_adjust(wspace=.25,hspace=.35)
        s = ('velocity_spherical_vela{0:d}_'
             'projection_{1:s}_galFrame_absCells.png'.format(galNum,stat))
        print(s)
        fig.savefig(s,dpi=300,bbox_inches='tight')

        plt.close(fig)


