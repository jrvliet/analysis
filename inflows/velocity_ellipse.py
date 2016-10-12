
'''
Plots the velocity ellipse
'''

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.linalg as sl
import scipy.stats as st


pd.options.mode.chained_assignment = None

def mkHist(ax,x,y,xlab,ylab,radial=False):
    stat= 'count'
    numbins = 50
    if radial==True:
        binrange = [[0,6],[-400,400]]
    else:
        binrange = [[-400,400],[-400,400]]

    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,x,
                    statistic=stat,bins=numbins,range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Counts',rotation=270,fontsize=12)


# File names
dataloc = '/home/jacob/research/velas/vela2b/vela27/'
dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotmat = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

# Filament resrictions
loT, hiT = 10**3.5, 10**4.5
zcut = 0
thetacut = 80
limits = pd.read_csv('cloudLimits.csv')

expns0 = range(20,50)
expns = [i/100. for i in expns0]

vrMin,vrMax = 1e6,0
vthetaMin,vthetaMax = 1e6,0
vphiMin,vphiMax = 1e6,0


for a in expns:
    
    print(a)

    # Read in Rvir
    with open(dataloc+rotmat.format(a)) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    # Read in the data
    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname, 'data')
    
    # Read in density cuto      
    nRange = limits[limits.expn==a]
    loN = 10**(nRange['loN'].values[0])
    hiN = 10**(nRange['hiN'].values[0])

    # Select out the gas
    tempInd = (df['temperature']>loT) & (df['temperature']<hiT) 
    densInd = (df['density']>loN) & (df['density']<hiN)
    spacInd = (df['zRot']>zcut) & (df['theta']<thetacut)
    cloud = df[tempInd & densInd & spacInd]

    cloud['rMod'] = cloud['r']/rvir

    # Calculate spherical velocities
    # These equations use theta is in the xy plane and phi is measured from z
    cloud['vr'] = (cloud['xRot']*cloud['vxRot'] + cloud['yRot']*cloud['vyRot'] +
                    cloud['zRot']*cloud['vzRot']) / cloud['r']
    cloud['thetadot'] = ( (cloud['xRot']*cloud['vyRot'] - cloud['yRot']*cloud['vxRot']) /
                         (cloud['xRot']**2 + cloud['yRot']**2 ) )
    cloud['vtheta'] = cloud['r']*np.sin(cloud['phi'])*cloud['thetadot']
    cloud['phidot'] = ( (-1 / np.sqrt( 1 - (cloud['zRot']/cloud['r'])**2)) *
                (cloud['vzRot']/cloud['r'] - cloud['zRot']*cloud['vr']/cloud['r']**2) )

    cloud['vphi'] = cloud['r']*cloud['phidot']

    if cloud['vr'].min()<vrMin:
        vrMin = cloud['vr'].min()
    if cloud['vr'].max()>vrMax:
        vrMax = cloud['vr'].max()
    if cloud['vphi'].min()<vphiMin:
        vphiMin = cloud['vphi'].min()
    if cloud['vphi'].max()>vphiMax:
        vphiMax = cloud['vphi'].max()
    if cloud['vtheta'].min()<vthetaMin:
        vthetaMin = cloud['vtheta'].min()
    if cloud['vtheta'].max()>vthetaMax:
        vthetaMax = cloud['vtheta'].max()

    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(15,10))
    # Reverse labels for phi and theta so the plot has the correct conventions
    mkHist(ax1,cloud['rMod'],cloud['vr'],'r','vr',True)
    mkHist(ax2,cloud['rMod'],cloud['vtheta'],'r','vphi',True)
    mkHist(ax3,cloud['rMod'],cloud['vphi'],'r','vtheta',True)
    mkHist(ax4,cloud['vr'],cloud['vtheta'],'vr','vphi')
    mkHist(ax5,cloud['vr'],cloud['vphi'],'vr','vtheta')
    mkHist(ax6,cloud['vtheta'],cloud['vphi'],'vphi','vtheta')



#    ax1.scatter(cloud['rMod'],cloud['vr'],marker='o',color='green',alpha=0.01)
#    ax2.scatter(cloud['rMod'],cloud['vtheta'],marker='o',color='green',alpha=0.01)
#    ax3.scatter(cloud['rMod'],cloud['vphi'],marker='o',color='green',alpha=0.01)
#    ax4.scatter(cloud['vr'],cloud['vtheta'],marker='o',color='green',alpha=0.01)
#    ax5.scatter(cloud['vr'],cloud['vphi'],marker='o',color='green',alpha=0.01)
#    ax6.scatter(cloud['vtheta'],cloud['vphi'],marker='o',color='green',alpha=0.01)
#    
#    for ax in [ax1,ax2,ax3]:
#        ax.set_xlabel('r [Rvir]')
#        ax.set_xlim([0,6])
#    
#    ax1.set_ylim([-400,400])
#    ax2.set_ylim([-400,400])
#    ax3.set_ylim([-400,400])
#
#    ax4.set_xlim([-400,400])
#    ax5.set_xlim([-400,400])
#    ax6.set_xlim([-400,400])
#    ax4.set_ylim([-400,400])
#    ax5.set_ylim([-400,400])
#    ax6.set_ylim([-400,400])
#
#    ax1.set_ylabel('v r [km/s]')
#    ax2.set_ylabel('v theta [km/s]')
#    ax3.set_ylabel('v phi [km/s]')
#
#    ax4.set_xlabel('v r [km/s]')
#    ax5.set_xlabel('v r [km/s]')
#    ax6.set_xlabel('v theta [km/s]')
#    ax4.set_ylabel('v theta [km/s]')
#    ax5.set_ylabel('v phi [km/s]')
#    ax6.set_ylabel('v phi [km/s]')
#
    redshift = 1./a - 1.
    ax2.set_title('a = {0:.3f} \t z = {1:.3f}'.format(a,redshift))
    
    fig.tight_layout()
    s = 'velocity_ellipse_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)

    plt.close(fig)

    
print('Vr limits:     {0:.3f} - {1:.3f}'.format(vrMin,vrMax))
print('Vtheta limits: {0:.3f} - {1:.3f}'.format(vthetaMin,vthetaMax))
print('Vphi limits:   {0:.3f} - {1:.3f}'.format(vphiMin,vphiMax))










