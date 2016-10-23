
'''
For each field plotted in density_cuts_analysis_simple,
this code fits a line to each density cut and outputs the 
parameters and the r^2 statistic
'''


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std





simples = [(-5.5,-5.0),(-5.0,-4.5),(-4.5,-4.0),
           (-4.0,-3.5),(-3.5,-3.0),(-3.0,-2.5)]

# Data names
loc = '/home/jacob/research/code/analysis/inflows/'
boxloc = '/home/jacob/research/velas/vela2b/vela27/'
densityFilename = 'density_cuts_parameters_{0:s}.h5'
tempLabels = ['cool','warm','hot']

# Read in data
groupsList = []
for tempLabel in tempLabels:

    filename = densityFilename.format(tempLabel)
    df = pd.read_hdf(loc+filename, 'data')

    # Cut out rows with low number of cells
    numCellLim = 1e4
    df.drop( df[df['numCells']<numCellLim].index, inplace=True)

    # Calculate combined standard deviations
    df['locStd'] = np.sqrt( (df['xRotStd']**2 + df['yRotStd']**2) / df['zRotStd']**2 )

    # Group by the density cuts
    groupsList.append( df.groupby(['loN','hiN']) )

header ='field temperature loN hiN slope intercept slopePvalue interceptPvalue rsquared'.split()
fields = ('speedStd stdDev valongStd vperpStd locStd snIIStd snIaStd snIImean snIamean' 
         ' rMean rStd nHmean nHStd tMean tStd rMeanMod speedMean valongMean vperpMean'
         ' numCells vStatMean vStatStd rMean vrStd vzRotMean vzRotStd vrhoRotMean vrhoRotStd').split()

data = np.zeros(len(header),dtype='a16')

for field in fields:

    print(field)
    for groups,tempLabel in zip(groupsList,tempLabels):
        print('\t',tempLabel)
        
        for key, group in groups:

            if key in simples:
                
                x = sm.add_constant(group['a'])
                y = group[field]
                model = sm.OLS(y,x)
                results = model.fit()
                
                p = np.zeros(len(header),dtype='a16')
                p[0] = field
                p[1] = tempLabel
                p[2] = key[0]
                p[3] = key[1]
                p[4] = results.params[1]
                p[5] = results.params[0]
                p[6] = results.pvalues[1]
                p[7] = results.pvalues[1]
                p[8] = results.rsquared

                data = np.vstack((data,p))

                
data = np.delete(data,(0),axis=0)
df = pd.DataFrame(data,columns=header)
outfile = 'density_cuts_lineFits.h5'
df.to_hdf(outfile,'data',mode='w')

                
                
            


