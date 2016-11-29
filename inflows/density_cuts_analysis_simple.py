
'''
A code to analyze the output of the density_cut code, a file
named projected_distance_distribution.h5
'''

from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mp



def mkPlot(ax,groups,fieldBase,stat,tempLabel):

    logfields = ['snIIStd','snIaStd','snIImean','snIamean',
             'nHmean','nHStd','nHmean','tMean','tStd',
             'vStatMean','vStatStd','numCells']
    logfields = 'SNII SNIa density temperature pressure'.split()
    simples = [(-5.5,-5.0),(-5.0,-4.5),(-4.5,-4.0),
           (-4.0,-3.5),(-3.5,-3.0),(-3.0,-2.5)]
    markers = ['o','s','^','*','x','D']

    for i, (densityKey,group) in enumerate(groups):
        if densityKey in simples:
            marker = markers[simples.index(densityKey)]
            ax.plot(group['a'],group[fieldBase+stat],label=densityKey,
                    linestyle='solid',marker=marker)
 

    if fieldBase in logfields and stat!='Ratio':
        ax.set_yscale('log')
    ax.set_xlabel('a')
    ax.set_ylabel('{0:s} {1:s}'.format(fieldBase,stat))
    ax.set_xlim([0.2,0.5])
    ax.set_title(tempLabel)

def mkLines(axes,lower,upper):

    #   First two are red and denote the start and end of the
    #   major merger
    #   Second two are green are denote the start and end of the
    #   "green zone", where some properties seem to diverge
    #   for currently unknown reasons
    for ax in axes:
        ax.vlines(0.29,lower,upper,colors='r',
                  linestyles='dashed',linewidth=1)
        ax.vlines(0.32,lower,upper,colors='r',
                  linestyles='dashed',linewidth=1)
        ax.vlines(0.35,lower,upper,colors='g',
                  linestyles='dashed',linewidth=1)
        ax.vlines(0.45,lower,upper,colors='g',
                  linestyles='dashed',linewidth=1)





excludeCGM = True

# Read in times
times = pd.read_csv('lookback_time.csv')
times = times.set_index('a')

# Read in data
loc = '/home/jacob/research/code/analysis/inflows/'
boxloc = '/home/jacob/research/velas/vela2b/vela27/'
boxname = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'

if excludeCGM:
    filebase = 'density_cuts_parameters_{0:s}_noCGM.h5'
else:
    filebase = 'density_cuts_parameters_{0:s}.h5'
tempLabels = ['cool','warm','hot']

groupsList = []
for tempLabel in tempLabels:

    filename = filebase.format(tempLabel)

    df = pd.read_hdf(loc+filename,'data')

    # Cut out rows with low number of cells
    numCellLim = 1e4
    df = df[df['numCells']>numCellLim]
    print('\n{0:s}'.format(tempLabel))
    print(df.columns)

    # Add a column to the dataframe that is the age of the
    # universe at that expnansion parameter
    age = np.zeros(len(df['a']))
    for i,a in enumerate(df['a']):
        index = float('{0:.2f}'.format(a))
        age[i] = times['age [Gyr]'].ix[index]
    df['age'] = age

    # Calculate combined standard deviations
    df['locStd'] = np.sqrt((df['xRotStd']**2 + df['yRotStd']**2)/df['zRotStd']**2)

    groups = df.groupby(['loN','hiN'])

    groupsList.append(groups)

print('Data read in')

# Properties to be plotted
fields = ['speedStd','stdDev','valongStd','vperpStd','locStd',
          'snIIStd','snIaStd','snIImean','snIamean',
          'rMean','rStd','nHmean','nHStd','tMean','tStd',
          'rMeanMod','speedMean','valongMean','vperpMean', 'numCells',
          'vStatMean','vStatStd','vrMean','vrStd',
          'vzRotMean','vzRotStd','vrhoRotMean','vrhoRotStd']
fields = ('speed valong vperp r rMod thetaRot phiRot '
          'SNII SNIa density temperature mass pressure '
          'vr vzRot vrhoRot vthetaRot vphiRot thermalV').split()
stats = 'Mean Std Ratio MeanMW Median'.split()

print(len(fields))
# Loop over the fields and plot them
for i in range(len(fields)):
    fieldBase = fields[i]
    print(fieldBase)

    for stat in stats:
        print('\t{0:s}'.format(stat))
        fig,axes = plt.subplots(1,3,figsize=(15,5))
        for groups,tempLabel,ax in zip(groupsList,tempLabels,axes):
            mkPlot(ax,groups,fieldBase,stat,tempLabel)

        lowest, uppest = axes[0].get_ylim()
        for ax in axes[1:]:
            lower, upper = ax.get_ylim()
            if lower<lowest:
                lowest = lower
            if upper>uppest:
                uppest = upper
        for ax in axes:
            ax.set_ylim([lowest,uppest])
        
        mkLines(axes,lowest,uppest)
        axes[0].legend(loc='upper left',ncol=1,labelspacing=0,
                frameon=True,fontsize='small')
        
        if excludeCGM:
            s = './denseCuts/simple/density_cuts_evolution_noCGM_{0:s}.png'.format(fieldBase+stat)
        else:
            s = './denseCuts/simple/density_cuts_evolution_{0:s}.png'.format(fieldBase+stat)
        fig.savefig(s,bbox_inches='tight',dpi=300)
        plt.close(fig)



