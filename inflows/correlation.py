
'''
Determiens correlations between various properties and
galactocentric distance
'''

from __future__ import print_function
import pandas as pd
import numpy as np

dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/a0.490/'
filename = 'vela2b-27_GZa0.490.rot.h5'


fname = dataloc+filename
df = pd.read_hdf(fname,'data')

df['r'] = np.sqrt(df['x']**2 + df['y']**2 + df['z']**2)
df['speed'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)

fields = 'speed density temperature SNII SNIa'.split()
logfields = 'density temperature SNII SNIa'.split()
for field in logfields:
    df['log'+field] = np.log10(df[field])

methods = 'pearson kendall spearman'.split()

cor = np.zeros((len(logfields)*2,len(methods)))


for i,field in enumerate(logfields):
    print('{0:s}'.format(field))
    for j,method in enumerate(methods):
        print('\t{0:s}'.format(method))
        cor[i,j] = df['r'].corr(df[field], method=method)

df = pd.DataFrame(cor,index=logfields,columns=methods)
outfile = 'spatialCorrelation_a0.490.h5'
df.to_hdf(outfile,'data')
        







