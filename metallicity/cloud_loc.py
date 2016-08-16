
from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np


baseloc = '/mnt/cluster/abs/cgm/vela2b/'
galnum = range(21,30)
expn = ['0.490']*len(galnum)
expn[galnum.index(24)] = '0.450'

fig1 = plt.figure(figsize=(15,15))
fig2 = plt.figure(figsize=(15,15))

filename =  '{0:s}vela{1:d}/a{2:s}/vela2b-{1:d}_GZa{2:s}.h5'
galcoords = '{0:s}vela{1:d}/a{2:s}/vela2b-{1:d}_a{2:s}_i90_cellCoords.h5'
i = 1

loN, hiN = 10**-5, 10**-4
loT, hiT = 10**4, 10**4.5

for galID, a in zip(galnum, expn):

    fname = filename.format(baseloc,galID,a)
    cname = galcoords.format(baseloc,galID,a)

    coords = pd.read_hdf(cname, 'data')
    box = pd.read_hdf(fname, 'data')

    xgal = coords['xgal']
    ygal = coords['ygal']
    zgal = coords['zgal']

    t = box['temperature']
    n = box['density']
    sn = box['SNII']
    x, y, z, metal = [], [], [], []
    for j in range(len(box['SNII'])):
        if t[j]<hiT and t[j]>loT and n[j]<hiN and n[j]>loN:
            x.append(xgal[j])
            y.append(ygal[j])
            z.append(zgal[j])
            metal.append(np.log10(sn[j]))

    print('Vela2b-{0:d}\n\t{1:d} cells in cloud'.format(galID, len(x)))

    ax = fig1.add_subplot(3, 3, i, projection='3d')
    #ax.scatter(d['x'], d['y'], d['z'], c='k', marker='.')
    sc = ax.scatter(x, y, z, c=metal, marker='o', alpha=0.1)
    ax.view_init(elev=90, azim=0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('vela2b-{0:d}'.format(galID))
    cbar = plt.colorbar(sc, ax=ax, use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Log MF SNII', rotation=270, fontsize=11)

    ax = fig2.add_subplot(3, 3, i, projection='3d')
    ax.scatter(x, y, z, c=metal, marker='o', alpha=0.1)
    ax.view_init(elev=0, azim=90)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('vela2b-{0:d}'.format(galID))
    cbar = plt.colorbar(sc, ax=ax, use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Log MF SNII', rotation=270, fontsize=11)

    i+=1

fig1.tight_layout()
fig1.savefig('vela2b_lowZ_loc_azim0_galcoords.png', bbox_inches='tight')

fig2.tight_layout()
fig2.savefig('vela2b_lowZ_loc_azim90_galcoords.png', bbox_inches='tight')


