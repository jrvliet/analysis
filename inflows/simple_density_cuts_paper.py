
'''
A code to analyze the output of the density_cut code, a file
named projected_distance_distribution.h5

For the paper!
'''

from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mp



def mkPlot(ax,groups,fieldBase,stat,tempLabel,ylab):

    
    logfields = 'SNII SNIa density temperature pressure'.split()
    markers = 'o ^ s'.split()
    colors = 'b g r'.split()
    labels = 'diffuse mid dense'.split()

    for i, (densityKey,group) in enumerate(groups):
        marker = markers[i]
        color = colors[i]
        redshift = 1./group['a']-1
        label = labels[i]
        print(densityKey,label)
        ax.plot(redshift,group[fieldBase+stat],label=label,
                linestyle='solid',marker=marker,color=color)


    if fieldBase in logfields and stat!='Ratio':
        ax.set_yscale('log')
    ax.invert_xaxis()
    ax.set_xlabel('Redshift')
    ax.set_ylabel(ylab)

def mkLines(ax):

    #   First two are red and denote the start and end of the
    #   major merger
    #   Second two are green are denote the start and end of the
    #   "green zone", where some properties seem to diverge
    #   for currently unknown reasons
    points = [2.57,2.13,1.86]
    ymin,ymax = ax.get_ylim()
    for point in points:
        ax.vlines(point,ymin,ymax,colors='k',
                  linestyles='dashed',linewidth=1)
    ax.set_ylim([ymin,ymax])




excludeCGM = False

# Read in times
times = pd.read_csv('lookback_time.csv')
times = times.set_index('a')

# Read in data
loc = '/home/jacob/research/code/analysis/inflows/'
boxloc = '/home/jacob/research/velas/vela2b/vela27/'
boxname = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
plotLoc = loc+'simpleDenseCuts/'

if excludeCGM:
    filebase = 'simple_density_cuts_parameters_{0:s}_noCGM.h5'
else:
    filebase = 'simple_density_cuts_parameters_{0:s}.h5'
tempLabels = ['cool','warm','hot']

# Read in data
tempLabel = 'cool'
filename = filebase.format(tempLabel)
df = pd.read_hdf(loc+filename,'data')

# Cut out rows with low number of cells
numCellLim = 1e4
print('\n{0:s}'.format(tempLabel))

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

print('Data read in')

# Properties to be plotted
fields = ['speedStd','stdDev','valongStd','vperpStd','locStd',
          'snIIStd','snIaStd','snIImean','snIamean',
          'rMean','rStd','nHmean','nHStd','tMean','tStd',
          'rMeanMod','speedMean','valongMean','vperpMean', 'numCells',
          'vStatMean','vStatStd','vrMean','vrStd',
          'vzRotMean','vzRotStd','vrhoRotMean','vrhoRotStd']
fields = ('speed valong vperp r rMod thetaRot phiRot numCellsFrac '
          'SNII SNIa density temperature mass pressure '
          'vr vzRot vrhoRot vthetaRot vphiRot thermalV').split()
stats = 'Mean Std Ratio MeanMW Median'.split()

fields = 'vr SNII'.split()
stats = 'Mean Std Median'.split()
ylabs = ['Radial Velocity [km/s]','SNII Mass Fraction']

print(len(fields))
# Loop over the fields and plot them
for i in range(len(fields)):
    fieldBase = fields[i]
    print(fieldBase)

    for stat in stats:
        fig,ax = plt.subplots(1,1,figsize=(5,5))
        mkPlot(ax,groups,fieldBase,stat,tempLabel,ylabs[i])

        #ax.set_ylim([lowest,uppest])
        
        mkLines(ax)
        ax.legend(loc='best',ncol=1,labelspacing=0,
                frameon=True)#,fontsize='small')
        if fieldBase=='SNII':
            ax.set_ylim([1e-5,1e-2])
        
        s = '{0:s}/simple_denseCut_{1:s}.pdf'.format(plotLoc,fieldBase+stat)
        fig.savefig(s,bbox_inches='tight',dpi=300)
        plt.close(fig)



