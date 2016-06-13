
'''
Calculates the EW distribution as a funciton of 
position angle of the LOS around the galaxy.
'''

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy import stats

def calc_pa(phi):
    '''
    Calculates the position angle for the LOS
    '''
    if phi<=90:
        return phi
    elif phi>90 and phi<=180:
        return 180.0 - phi
    elif phi>180 and phi<=270:
        return phi - 180.0
    else:
        return 350.0 - phi

def gmean(vals):
    '''
    Returns geometric mean of values in vals
    '''
    return stats.gmean(vals)

coplaneCut = np.radians(20)
outflowCut = np.radians(60)
numbins = 50

location = '/mnt/cluster/abs/cgm/vela2b/'
halo_list = [25, 26, 27, 28]
halo_list = [28]
expn = '0.490'
incline = 90
ions = ['HI', 'MgII', 'CIV', 'OVI']

absfileBase = 'vela2b-{0:d}.{1:s}.{2:s}.i{3:d}.abs_cells.h5'
sysabsfileBase = 'vela2b-{0:d}.{1:s}.a{2:s}.ALL.sysabs.h5'
stat = 'mean' 
statname = 'mean'

# Loop over galaxies

for galNum in halo_list:

    fig, ((ax1, ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(10,10))
    axes = (ax1, ax2, ax3, ax4)

    for ion, ax in zip(ions, axes):
        
        loc = '{0:s}/vela{1:d}/a{2:s}/i{3:d}/{4:s}/'.format(location,galNum,expn,incline,ion)
        absfile = absfileBase.format(galNum,expn,ion,incline)
        sysabsfile = sysabsfileBase.format(galNum,ion,expn)

        dcell = pd.read_hdf(loc+absfile, 'data')
        dabs = pd.read_hdf(loc+sysabsfile, 'data')

        phi = []
        for p in dabs['phi']:
            phi.append(calc_pa(p))


        ew = dabs['EW_r'][ dabs['EW_r']>0 ]
        logN = dabs['logN'][ dabs['logN']>0 ]

        b = dabs['D']

        x = [ i*np.cos(j) for i,j in zip(b,phi)]
        y = [ i*np.sin(j) for i,j in zip(b,phi)]
    
        # Make 2d histogram of x, y position on sky, with bin value
        bin_means, x_edges, y_edges, binnumber = stats.binned_statistic_2d(x,y,dabs['EW_r'],
                statistic='mean', bins=numbins)
        whereNans = np.isnan(bin_means)
        bin_means[whereNans] = 0.0
        bin_means = np.ma.masked_where(bin_means==0, bin_means)
        bin_means = np.log10(bin_means)
        mesh = ax.pcolormesh(x_edges, y_edges, bin_means)
        ax.set_xlim([0,x_edges[-1]])
        ax.set_ylim([0,y_edges[-1]])
        ax.set_title(ion)
        cbar = fig.colorbar(mesh, ax=ax, use_gridspec=True)
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel('Mean EW', rotation=270, fontsize=12)

        # Draw lines seperating out the coplanar regions and the outflow regions
        x = [0, x_edges[-1]*np.cos(coplaneCut)]
        y = [0, x_edges[-1]*np.sin(coplaneCut)]
        ax.plot(x,y,linestyle='dotted',color='k')
        x = [0, x_edges[-1]*np.cos(outflowCut)]
        y = [0, x_edges[-1]*np.sin(outflowCut)]
        ax.plot(x,y,linestyle='dotted',color='k')



    fig.tight_layout()
    fig.savefig('pahist_vela2b-{0:d}_ew_{1:s}.png'.format(galNum, statname), bbox_inches='tight')





