
'''
Calculates SFR for vela2b-27
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



dataloc = '/mnt/cluster/abs/Simulations/vela2.1/VELA27/output/ana/'
filename = 'halos_{0:.3f}.txt'

expns = np.arange(0.200,0.500,0.01)

for a in expns:

    fname = dataloc + filename.format(a)
    with open(fname) as f:

        f.readline()
        f.readline()
        l = f.readline().split()
        mstar = float(l[13])
        mvir = float(l[7])
        rvir = float(l[8])
        mgas = float(l[12])

    



