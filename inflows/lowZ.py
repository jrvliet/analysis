
'''
Selects out gas with lowest metallicity and plots location,
prints properties
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as st
pd.options.mode.chained_assignment = None


def mkHist(ax,x,y,xlabel,ylabel):
    numbins = 50
    stat = 'count'
    if 'z' in ylabel:
        binrange = [[-3,3],[0,6]]
    else:
        binrange = [[-3,3],[-3,3]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,y,
                              statistic=stat, bins=numbins,
                                 range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)

    mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)


    return h,xedges,yedges,binnumber



dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expn = np.arange(0.200,0.500,0.01)

header = 'a xRot yRot zRot r theta phi density temperature SNII SNIa'.split()
logfields = 'density temperature SNII SNIa'.split()

data = np.zeros((len(expn),len(header)))
datadf = pd.DataFrame(data,columns=header)

for i,a in enumerate(expn):

    print(a)
    redshift = 1./a - 1

    rname = dataloc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname,'data')

    # Get the value of SNII that cuts out bottom 10%
    snIILim = df['SNII'].quantile(q=0.1)
    metalLims = (df['SNII']<=snIILim)
#    lowZ = df[metalLims]
    spaceLims = (df['theta']<80)
    lowZ = df[metalLims & spaceLims]

    fig,axes = plt.subplots(1,3,figsize=(15,5))
    mkHist(axes[0],lowZ['xRot']/rvir,lowZ['yRot']/rvir,'xRot','yRot')
    mkHist(axes[1],lowZ['xRot']/rvir,lowZ['zRot']/rvir,'xRot','zRot')
    mkHist(axes[2],lowZ['yRot']/rvir,lowZ['zRot']/rvir,'yRot','zRot')
    axes[1].set_title('a={0:.3f}, z={1:.3f}'.format(a,redshift))

    fig.tight_layout()
    s = 'lowZ_spatial_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)
    
    datadf.iloc[i]['a'] = a
    for j,field in enumerate(header[1:]):
        if field in logfields:
            datadf.iloc[i][field] = np.log10(lowZ[field].mean())
        else:
            datadf.iloc[i][field] = lowZ[field].mean()
    
    
#datadf.to_hdf('lowZProps.h5','data',mode='w')
datadf.to_hdf('lowZProps_filament.h5','data',mode='w')




