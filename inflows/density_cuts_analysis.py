
'''
A code to analyze the output of the density_cut code, a file
named projected_distance_distribution.h5
'''

from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
import matplotlib as mp

print(mp.matplotlib_fname())

# Read in times
times = pd.read_csv('lookback_time.csv')
times = times.set_index('a')

# Read in data
loc = '/home/jacob/research/code/analysis/inflows/'
boxloc = '/home/jacob/research/velas/vela2b/vela27/'
boxname = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'

filebase = 'density_cuts_parameters_{0:s}.h5'
tempLabels = ['cool','warm','hot']


# Properties to be plotted
fields = ('speedStd stdDev valongStd vperpStd locStd '
          'snIIStd snIaStd snIImean snIamean '
          'rMean rStd nHmean nHStd tMean tStd '
          'rMeanMod speedMean valongMean vperpMean numCells '
          'vrMean vrStd vzRotMean vzRotStd vrhoRotMean vrhoRotStd '

fields = ('speed valong vperp r rMod thetaRot phiRot '
          'SNII SNIa density temperature mass pressure '
          'vr vzRot vrhoRot vthetaRot vphiRot thermalV').split()   


# Plotting options
logfields = ['snIIStd','snIaStd','snIImean','snIamean',
             'nHmean','nHStd','nHmean','tMean','tStd',
             'vStatMean','vStatStd']
logfields = 'SNII SNIa density temperature'.split()
upperlims = [150, 60, 150, 200, 4.00, 
             1e-1, 1e-2, 1e-1, 1e-2, 
             500, 250 , 1e-1, 1e-1, 1e5, 10**4.5,
             6, 600, 250, 500, 700000,
             1e4, 1e6, 600, 300,
             600, 600, 600, 600]
lowerlims = [0, 0, 0, 0, 0.10, 
             1e-5, 1e-8, 1e-5, 1e-8, 
             1, 1, 1e-6, 1e-6, 1e4, 1e3,
             0, 30, 0, 10, 0,
             0.1, 0.1, -300, 0, 
             -600, -600, -600, -600]
lines = ['-','--','-.']
markers = ['o','s','^','*','x']


for tempLabel in tempLabels:

    filename = filebase.format(tempLabel)
    df = pd.read_hdf(loc+filename,'data')

    # Cut out rows with low number of cells
    numCellLim = 1e4
    df = df[df['numCells']>numCellLim]

    # Add a column to the dataframe that is the age of the
    # universe at that expnansion parameter
    age = np.zeros(len(df['a']))
    for i,a in enumerate(df['a']):
        index = float('{0:.2f}'.format(a))
        age[i] = times['age [Gyr]'].ix[index]
    df['age'] = age

    # Calculate combined standard deviations
    df['locStd'] = np.sqrt((df['xRotStd']**2 + df['yRotStd']**2)/df['zRotStd']**2)
    #df['totStd'] = np.sqrt(df['locStd']**2 + df['speedStd']**2)

    print(np.log10(df['vStatMean'].min()))
    print(np.log10(df['vStatStd'].min()))

    groups = df.groupby(['loN','hiN'])

    simples = [(-5.5,-5.0),(-5.0,-4.5),(-4.5,-4.0),
               (-4.0,-3.5),(-3.5,-3.0),(-3.0,-2.5)]

    # Loop over various cuts of the number of cells 
    cellLims = np.arange(4.0,6,0.5)
    cellLims = [4.0]
    #dfFull = df.copy()
    for numCellLim in cellLims:
        print(numCellLim)

        # Select out cells with at least numCellLim cells
        df = df[df['numCells']>10**numCellLim]


        # Group the data set by combinations of low and high
        # density cuts
        groups = df.groupby(['loN','hiN'])
        
        # Loop over the fields and plot them
        for i in range(len(fields)):
            field = fields[i]

            figAll,axAll = plt.subplots(1,1,figsize=(10,10))
            figSimp,axSimp = plt.subplots(1,1,figsize=(10,10))
            figAll2,axAll2 = plt.subplots(1,1,figsize=(10,10))
            figSimp2,axSimp2 = plt.subplots(1,1,figsize=(10,10))

            linecycler = cycle(lines)
            markercycler = cycle(markers)
            for key, group in groups:
                axAll.plot(group['a'],group[field],label=key,
                        linestyle=next(linecycler),
                        marker=next(markercycler))
                axAll2.plot(group['age'],group[field],label=key,
                        linestyle=next(linecycler),
                        marker=next(markercycler))
                if key in simples:
                    axSimp.plot(group['a'],group[field],
                            label=key,
                            linestyle='solid',
                            marker=next(markercycler))
                    axSimp2.plot(group['age'],group[field],
                            label=key,
                            linestyle='solid',
                            marker=next(markercycler))
            
            for ax in [axAll,axSimp]:
                ax.set_xlabel('a')
                ax.set_ylabel(field)
                ax.set_ylim([lowerlims[i],upperlims[i]])
                ax.set_xlim([0.2,0.5])

                # Plot four vertical lines. 
                #   First two are red and denote the start and end of the
                #   major merger
                #   Second two are green are denote the start and end of the
                #   "green zone", where some properties seem to diverge
                #   for currently unknown reasons
                ax.vlines(0.29,lowerlims[i],upperlims[i],colors='r',
                          linestyles='dashed',linewidth=2)
                ax.vlines(0.32,lowerlims[i],upperlims[i],colors='r',
                          linestyles='dashed',linewidth=2)
                ax.vlines(0.35,lowerlims[i],upperlims[i],colors='g',
                          linestyles='dashed',linewidth=2)
                ax.vlines(0.45,lowerlims[i],upperlims[i],colors='g',
                          linestyles='dashed',linewidth=2)
                if field in logfields:
                    ax.set_yscale('log')
            for ax in [axAll2,axSimp2]:
                ax.set_xlabel('Age [Gyr]')
                ax.set_ylabel(field)
                ax.set_ylim([lowerlims[i],upperlims[i]])
                ax.set_xlim([1.5,6.0])

                # Plot four vertical lines. 
                #   First two are red and denote the start and end of the
                #   major merger
                #   Second two are green are denote the start and end of the
                #   "green zone", where some properties seem to diverge
                #   for currently unknown reasons
                ax.vlines(2.68,lowerlims[i],upperlims[i],colors='r',
                          linestyles='dashed',linewidth=2)
                ax.vlines(3.10,lowerlims[i],upperlims[i],colors='r',
                          linestyles='dashed',linewidth=2)
                ax.vlines(3.54,lowerlims[i],upperlims[i],colors='g',
                          linestyles='dashed',linewidth=2)
                ax.vlines(5.07,lowerlims[i],upperlims[i],colors='g',
                          linestyles='dashed',linewidth=2)
                if field in logfields:
                    ax.set_yscale('log')

            axAll.legend(loc='upper left',ncol=5,labelspacing=0,frameon=True)
            axSimp.legend(loc='upper left',ncol=1,labelspacing=0,frameon=True)
            axAll2.legend(loc='upper left',ncol=5,labelspacing=0,frameon=True)
            axSimp2.legend(loc='upper left',ncol=1,labelspacing=0,frameon=True)
            s = './denseCuts/density_cuts_evolution_{0:s}_{1:.1f}_{2:s}.png'.format(field,numCellLim,tempLabel)
            figAll.savefig(s,bbox_inches='tight',dpi=300)
            s = './denseCuts/density_cuts_evolution_{0:s}_{1:.1f}_{2:s}_time.png'.format(field,numCellLim,tempLabel)
            figAll2.savefig(s,bbox_inches='tight',dpi=300)
            s = './denseCuts/simple_density_cuts_evolution_{0:s}_{1:.1f}_{2:s}.png'.format(field,numCellLim,tempLabel)
            figSimp.savefig(s,bbox_inches='tight',dpi=300)
            s = './denseCuts/simple_density_cuts_evolution_{0:s}_{1:.1f}_{2:s}_time.png'.format(field,numCellLim,tempLabel)
            figSimp2.savefig(s,bbox_inches='tight',dpi=300)

            plt.close(figAll)
            plt.close(figSimp)
            plt.close(figAll2)
            plt.close(figSimp2)



