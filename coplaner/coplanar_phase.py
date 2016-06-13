
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
    print h.max()

    h = np.rot90(h)
    h = np.flipud(h)
    h = np.ma.masked_where(h==0, h)

    h = np.log10(h)

    print h.max()
    # Generate the plots
    cmap = plt.cm.get_cmap('viridis')
    mesh = ax.pcolormesh(xedges, yedges, h, cmap=cmap)
    ax.set_xlim([-8,1])
    ax.set_ylim([2,8])
    ax.set_xlabel('$\log (n_{H})$')
    ax.set_ylabel('$\log (T)$')
    ax.set_title(label, fontsize=10)
    
    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    #cbar.ax.set_ylabel(cbarLabel, rotation=270, fontsize=12)


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
perpCut = 60       # PA below which the cells are considered coplanar
impactCut = 0.1    # Impact papameter below which it will be considered ISM [Rvir]
numbins = 50
binrange = [[-8,1],[2,8]]

ions = ['MgII','CIV', 'OVI']

numions = len(ions)
fig, axes = plt.subplots(numions, 3, figsize=(10,10))

for ion, axrow in zip(ions, axes):
    testloc = '/mnt/cluster/abs/cgm/vela2b/vela28/a0.490/i90/{0:s}/'.format(ion)
    testfile = 'vela2b-28.0.490.{0:s}.i90.abs_cells.h5'.format(ion)
    filename = '{0:s}/{1:s}'.format(testloc,testfile)

    # Get the virial radius
    with open(testloc+'../galaxy.props', 'r') as f:
        for i in range(5):
            line = f.readline()
        rvir = float(line.split()[1])

    # Read in the file
    d = pd.read_hdf(filename, 'data')
    numcells = d.shape[0]

    # Read in lines information
    losnum, b, phi = np.loadtxt(testloc+'lines.info', skiprows=1, usecols=(0,1,2), unpack=True)

    # Convert the angle to position angle
    pa = convert_angle(phi)

    # Coplanar density and temperature
    coDense, coTemp = [], []
    nonDense, nonTemp = [], []
    perpDense, perpTemp = [], []

    # Loop over the absorption cells
    for i in range(1, numcells):
        
        index = int(d['LOS'][i]) - 1
        angle = pa[index]
        impact = b[index]
        
        if impact > 0.1*rvir:
            if angle < coplanarCut:
                coDense.append( d['log_nH'][i] )
                coTemp.append( d['log_T'][i] )
            elif angle > perpCut:
                perpDense.append( d['log_nH'][i] )
                perpTemp.append( d['log_T'][i] )
            else:
                nonDense.append( d['log_nH'][i] )
                nonTemp.append( d['log_T'][i] )
            
    # Bin the data and plot
    label = '{0:s} Within {1:d} deg of Major Axis'.format(ion, coplanarCut)
    plot_hist(coDense, coTemp, numbins, binrange, axrow[0], label)

    label = '{0:s} In Void'.format(ion)
    plot_hist(nonDense, nonTemp, numbins, binrange, axrow[1], label)

    label = '{0:s} Within {1:d} deg of Minor Axis'.format(ion, 90-perpCut)
    plot_hist(perpDense, perpTemp, numbins, binrange, axrow[2], label)


fig.tight_layout()
s = 'vela2b-28_a0.490_i90_EWcut100mA.pdf'
plt.savefig(s, bbox_inches='tight')




