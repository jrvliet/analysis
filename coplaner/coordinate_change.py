
import pandas as pd
import numpy as np
import os



galIDs = range(21,30)
expns = ['0.490']*len(galIDs)
expns[galIDs.index(24)] = '0.450'

for galID, a in zip(galIDs, expns):
    print galID

    loc = './vela{0:d}/a{1:s}/'.format(galID,a)
    gasbox = 'vela2b-{0:d}_GZa{1:s}.h5'.format(galID, a)
    d = pd.read_hdf(loc+gasbox, 'data')

    # Read in the rotation matrix
    rotmat = np.zeros((3,3))
    rotmatFile = loc+'rotmat_a{0:s}.txt'.format(a)
    with open(rotmatFile, 'r') as f:
        f.readline()
        l = f.readline().split()
        rvir = float(l[3])
        a11 = float(l[4])
        a12 = float(l[5])
        a13 = float(l[6])
        a21 = float(l[7])
        a22 = float(l[8])
        a23 = float(l[9])
        a31 = float(l[10])
        a32 = float(l[11])
        a33 = float(l[12])

    galcoordFile = 'vela2b-{0:d}_a{1:s}_galaxy.coords'.format(galID,a)
    isfile = os.path.isfile(galcoordFile)
    isfile = 0
    if not isfile:
        s = '{0:f}\t{1:f}\t{2:f}\n'
        with open(galcoordFile, 'w') as f:
            for i in range(len(d['x'])):
                x = a11*d['x'][i] + a12*d['y'][i] + a13*d['z'][i]
                y = a21*d['x'][i] + a22*d['y'][i] + a23*d['z'][i]
                z = a31*d['x'][i] + a32*d['y'][i] + a33*d['z'][i]
                f.write(s.format(x, y, z))

