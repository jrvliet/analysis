
'''
Determiens correlations between various properties and
galactocentric distance
'''

import pandas as pd
import numpy as np

dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/a0.490/'
filename = 'vela2b-27_GZa0.490.rot.h5'


fname = dataloc+filename
df = pd.read_hdf(fname,'data')

df['r'] = np.loadtxt(df['x']**2 + df['y']**2 + df['z']**2)

logfields = 'density temperature SNII SNIa'.split()
for field in logfields:
    df[field] = np.log10(df[field])

methods = 'pearson kendall spearman'.split()

cor = np.zeros((len(logfields),len(methods)))


for i,field in enumerate(logfields):
    for j,method in enumerate(methods):

        cor[i,j] = df['r'].corr(field, method=method)

df = pd.DataFrame(cor,index=logfields,columns=method)
outfile = 'spatialCorrelation_a0.490.h5'
df.to_hdf(outfile,'data')
        







