
'''
Calculates the phase of all absorbing cells, seperated into
cells that are coplanar and those that are not
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys



def plot_hist(x, y, numbins, binrange, ax, label):

    h, xedges, yedges = np.histogram2d( x, y, bins=numbins, range=binrange)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    h = np.rot90(h)
    h = np.flipud(h)
    h = np.ma.masked_where(h==0, h)

    print label, np.sum(h)
    h = np.log10(h)

    # Generate the plots
    cmap = plt.cm.get_cmap('viridis')
    mesh = ax.pcolormesh(xedges, yedges, h, cmap=cmap)
    ax.set_xlim([-8,1])
    ax.set_ylim([2,8])
    ax.set_xlabel('$\log (n_{H})$')
    ax.set_ylabel('$\log (T)$')
    ax.set_title(label, fontsize=10)
    
    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbarLabel = 'Counts'
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(cbarLabel, rotation=270, fontsize=12)


def convert_angle(phi):

    pa = []
    
    for angle in phi:

        if angle<=90:
            pa.append(angle)
        elif angle>90 and angle<=180:
            pa.append(180-angle)
        elif angle>180 and angle<=270:
            pa.append(angle-180)
        else:
            pa.append(360-angle)
    
    return pa




coplanarCut = 20   # PA below which the cells are considered coplanar
perpCut = 45       # PA above which the cells are considered outflow
impactCut = 0.2    # Impact papameter below which it will be considered ISM [Rvir]
ewCut = 0.1        # EW cut [A]

planeSize = coplanarCut
voidSize = perpCut - coplanarCut
minorSize = 90 - perpCut

numbins = 50
binrange = [[-8,1],[2,8]]

ions = ['MgII','CIV', 'OVI']
galID = 28
expn = '0.490'
inc = 90

numions = len(ions)
fig, axes = plt.subplots(numions, 4, figsize=(12,9))

for ion, axrow in zip(ions, axes):
    testloc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:s}/i{2:d}/{3:s}/'.format(galID, expn, inc,ion)
    absfile = 'vela2b-{0:d}.{1:s}.{2:s}.i{3:d}.abs_cells.h5'.format(galID, expn, ion, inc)
    sysfile = 'vela2b-{0:d}.{1:s}.a{2:s}.i{3:d}.ALL.sysabs.h5'.format(galID,ion,expn,inc)
    filename = '{0:s}/{1:s}'.format(testloc,absfile)

    # Get the virial radius
    with open(testloc+'../galaxy.props', 'r') as f:
        for i in range(5):
            line = f.readline()
        rvir = float(line.split()[1])

    # Read in the files
    d = pd.read_hdf(filename, 'data')
    numcells = d.shape[0]
    dsys = pd.read_hdf(testloc+sysfile, 'data')
    ew = dsys['EW_r']

    # Read in lines information
    losnum, b, phi = np.loadtxt(testloc+'lines.info', skiprows=1, usecols=(0,1,2), unpack=True)

    # Convert the angle to position angle
    pa = convert_angle(phi)

    # Coplanar density and temperature
    coDense, coTemp = [], []
    nonDense, nonTemp = [], []
    perpDense, perpTemp = [], []
    dense, temp = [], []

    # Loop over the absorption cells
    for i in range(1, numcells):
        
        index = int(d['LOS'][i]) - 1
        angle = pa[index]
        impact = b[index]
        width = ew[index]
        
        if impact > impactCut*rvir and width>=ewCut:
            if angle < coplanarCut:
                coDense.append( d['log_nH'][i] )
                coTemp.append( d['log_T'][i] )
            elif angle > perpCut:
                perpDense.append( d['log_nH'][i] )
                perpTemp.append( d['log_T'][i] )
            else:
                nonDense.append( d['log_nH'][i] )
                nonTemp.append( d['log_T'][i] )
            dense.append( d['log_nH'][i])
            temp.append( d['log_T'][i])
            
    # Bin the data and plot
    label = '{0:s} Within {1:d} deg of Major Axis'.format(ion, coplanarCut)
    plot_hist(coDense, coTemp, numbins, binrange, axrow[0], label)

    label = '{0:s} In Void'.format(ion)
    plot_hist(nonDense, nonTemp, numbins, binrange, axrow[1], label)

    label = '{0:s} Within {1:d} deg of Minor Axis'.format(ion, 90-perpCut)
    plot_hist(perpDense, perpTemp, numbins, binrange, axrow[2], label)

    label = '{0:s} All'.format(ion)
    plot_hist(dense, temp, numbins, binrange, axrow[3], label)


fig.tight_layout()
s = 'vela2b-{0:d}_a{1:s}_i{2:d}_EWcut{3:d}mA_PAphase.pdf'.format(galID,expn,inc,int(ewCut*1000))
plt.savefig(s, bbox_inches='tight')




