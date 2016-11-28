

''' 
Determines the escape velocity for finely meshed grid of distances from the host
galaxy
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt



dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expn = np.arange(0.200,0.500,0.01)
rmin = 0
rmax = 5
numRbins = 50
rbins = np.
rbins = np.linspace(rmin,rmax,numRbins+1)

escapeProper = np.zeros((len(rbins),len(expn)))
escapeRvir = np.zeros((len(rbins),len(expn)))

for a in expn:

    rname = dataloc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname,'data')

    df['r'] = np.sqrt(df['x']**2 + df['y']**2 + df['z']**2)


    for 




# This won't work, need stars and DM as well -> look for mass profiles





