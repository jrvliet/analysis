
from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



baseloc = '/mnt/cluster/abs/cgm/vela2b/'
galnum = range(21,30)
expn = ['0.490']*len(galnum)
expn[galnum.index(24)] = '0.450'


for galID, a in zip(galnum, expn):

    loc = '{0:s}vela{1:d}/a{2:s}/'.format(baseloc, galID, a)
    filename = '{0:s}vela2b-{1:d}_GZa{2:s}.h5'.format(loc,galID,a)

    d = pd.read_hdf(filename, 'data')

    lowCut = d.quantile(0.1)['SNII']

    losZd = d[(d.SNII<lowCut)]

    r = d['x']**2



