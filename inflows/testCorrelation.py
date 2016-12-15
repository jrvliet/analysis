
from __future__ import print_function
import numpy as np
import pandas as pd

loc = '/home/jacob/research/velas/vela2b/vela27/'
loc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
expns = [0.540]
expns = np.arange(0.200,0.550,0.01)

hdfFile = 'distanceCorrelation.h5'
store = pd.HDFStore(hdfFile)

fields = 'r speed density temperature SNII SNIa'.split()
logfields = 'logdensity logtemperature logSNII logSNIa'.split()
methods = 'pearson kendall spearman'.split()
methods = 'pearson spearman'.split()

allfields = fields + logfields

for a in expns:
    print(a)

    cor = np.zeros((len(allfields),len(methods)))
    cordf = pd.DataFrame(cor,columns=methods,index=allfields)

    fname = loc+filename.format(a)
    
    df = pd.read_hdf(fname,'data')
    
    df['r'] = np.sqrt(df['x']**2 + df['y']**2 + df['z']**2)
    df['speed'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)

    for logfield in logfields:
        basefield = logfield.split('log')[-1]
        df[logfield] = np.log10(df[basefield])

    d = df[allfields]

    for method in methods:
        print('\tMethod = {0:s}'.format(method))
        cordf[method] = d.corr(method=method)['r']


    storeName = 'a{0:.3f}'.format(a)
    storeName = 'a{0:d}'.format(int(a*1000))
    store[storeName] = cordf

store.close()
