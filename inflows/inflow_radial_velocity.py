
'''
Plots the radial velocity of cells in the inflows
of vela2b-27 as a funciton of galactocentric distance
'''

from __future__ import print_function
import pandas as np
import numpy as np
import matplotlib.pyplot as plt

fname = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'

expns = np.arange(0.200,0.510,0.01)

for a in expns:

    fn = loc.format(a) + fname.format(a)

    










