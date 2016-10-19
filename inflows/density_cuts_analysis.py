
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
filename = 'projected_distance_distribution.h5'
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

# Properties to be plotted
fields = ['speedStd','stdDev','valongStd','vperpStd','locStd',
          'snIIStd','snIaStd','snIImean','snIamean',
          'rMean','rStd','nHmean','nHStd','tMean','tStd',
          'rMeanMod','speedMean','valongMean','vperpMean', 'numCells',
          'vStatMean','vStatStd','vrMean','vrStd']

# Plotting options
logfields = ['snIIStd','snIaStd','snIImean','snIamean',
             'nHmean','nHStd','nHmean','tMean','tStd',
             'vStatMean','vStatStd']
upperlims = [80, 60, 90, 70, 4.00, 
             1e-1, 1e-2, 1e-1, 1e-2, 
             500, 250 , 1e-1, 1e-1, 1e5, 10**4.5,
             6, 300, 250, 175, 700000,
             1e4, 1e6, 200, 200]
lowerlims = [0, 0, 0, 0, 0.10, 
             1e-5, 1e-8, 1e-5, 1e-8, 
             1, 1, 1e-6, 1e-6, 1e4, 1e3,
             0, 30, 0, 10, 0,
             0.1, 0.1, -300, 0]
lines = ['-','--','-.']
markers = ['o','s','^','*','x']
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
        s = './denseCuts/density_cuts_evolution_{0:s}_{1:.1f}.png'.format(field,numCellLim)
        figAll.savefig(s,bbox_inches='tight',dpi=300)
        s = './denseCuts/density_cuts_evolution_{0:s}_{1:.1f}_time.png'.format(field,numCellLim)
        figAll2.savefig(s,bbox_inches='tight',dpi=300)
        s = './denseCuts/simple_density_cuts_evolution_{0:s}_{1:.1f}.png'.format(field,numCellLim)
        figSimp.savefig(s,bbox_inches='tight',dpi=300)
        s = './denseCuts/simple_density_cuts_evolution_{0:s}_{1:.1f}_time.png'.format(field,numCellLim)
        figSimp2.savefig(s,bbox_inches='tight',dpi=300)

        plt.close(figAll)
        plt.close(figSimp)
        plt.close(figAll2)
        plt.close(figSimp2)



## In[ ]:
#
#loloc, hiloc, aloc = [], [], []
#lokin, hikin, akin = [], [], []
#lotot, hitot, atot = [], [], []
#for key,group in groups:
#    # Location based
#    minIndex = group['locStd'].idxmin()
#    loN = group['loN'].ix[minIndex]
#    hiN = group['hiN'].ix[minIndex]
#    loloc.append(loN)
#    hiloc.append(hiN)
#    aloc.append(key)
#    # Kinematics based
#    minIndex = group['speedStd'].idxmin()
#    loN = group['loN'].ix[minIndex]
#    hiN = group['hiN'].ix[minIndex]
#    lokin.append(loN)
#    hikin.append(hiN)
#    akin.append(key)
#    # Both based
#    minIndex = group['totStd'].idxmin()
#    loN = group['loN'].ix[minIndex]
#    hiN = group['hiN'].ix[minIndex]
#    lotot.append(loN)
#    hitot.append(hiN)
#    atot.append(key)
#endloclo = loloc[-1]
#endlochi = hiloc[-1]
#endtotlo = lotot[-1]
#endtothi = hitot[-1]
#
#
## In[ ]:
#
#fig,ax = plt.subplots(1,1,figsize=(10,10))
#ax.plot(aloc,loloc,color='r',label='Loc')
#ax.plot(aloc,hiloc,color='r')
#
#ax.plot(akin,lokin,color='b',label='Kin')
#ax.plot(akin,hikin,color='b')
#
#ax.plot(atot,lotot,color='g',label='Tot')
#ax.plot(atot,hitot,color='g')
#
#ax.legend(loc='upper left')
#ax.set_ylim([-6,-2])
#
#
## In[ ]:
#
#keys = [0.49]
#for key in keys:
#    rotmat = boxloc+'a{0:.3f}/rotmat_a{0:.3f}.txt'.format(key)
#    with open(rotmat) as f:
#        f.readline()
#        rvir = float(f.readline().split()[3])
#    fname = boxloc+boxname.format(key)
#    df = pd.read_hdf(fname,'data')
#    loT, hiT = 10**3.5, 10**4.5
#    lowN, highN = endloclo, endlochi
#    loN = 10**lowN
#    hiN = 10**highN
#    tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
#    densInd = (df['density']>loN) & (df['density']<hiN)
#    spacInd = (df['zRot']>0) & (df['theta']<80)
#    cloud = df[tempInd & densInd & spacInd]
#    mkPlot(cloud,rvir,key,lowN,highN,'loc')
#    
#    df = pd.read_hdf(fname,'data')
#    loT, hiT = 10**3.5, 10**4.5
#    lowN, highN = endtotlo, endtothi
#    loN = 10**lowN
#    hiN = 10**highN
#    tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
#    densInd = (df['density']>loN) & (df['density']<hiN)
#    spacInd = (df['zRot']>0) & (df['theta']<80)
#    cloud = df[tempInd & densInd & spacInd]
#    mkPlot(cloud,rvir,key,lowN,highN,'tot')
#
#
## In[ ]:
#
#boxloc = '/home/jacob/research/velas/vela2b/vela27/'
#boxname = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
#for key, group in groups:
#    print(key)
#    rotmat = boxloc+'a{0:.3f}/rotmat_a{0:.3f}.txt'.format(key)
#    with open(rotmat) as f:
#        f.readline()
#        rvir = float(f.readline().split()[3])
#    fname = boxloc+boxname.format(key)
#    df = pd.read_hdf(fname,'data')
#    loT, hiT = 10**3.5, 10**4.5
#    #minIndex = group['ratio'].idxmin()
#    #lowN = group['loN'].ix[minIndex]
#    #highN = group['hiN'].ix[minIndex]
#    lowN = -5.5
#    highN = -3.5
#    loN = 10**lowN
#    hiN = 10**highN
#    tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
#    densInd = (df['density']>loN) & (df['density']<hiN)
#    spacInd = (df['zRot']>0) & (df['theta']<80)
#    cloud = df[tempInd & densInd & spacInd]
#    mkPlot(cloud,rvir,key,lowN,highN)
#
#
## In[ ]:
#
#def mkPlot(cl,rvir,a,loN,hiN,label=''):
#    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
#    ax1.scatter(cl['xRot']/rvir,cl['yRot']/rvir,marker='o',color='green',alpha=0.01)
#    ax2.scatter(cl['xRot']/rvir,cl['zRot']/rvir,marker='o',color='green',alpha=0.01)
#    ax3.scatter(cl['yRot']/rvir,cl['zRot']/rvir,marker='o',color='green',alpha=0.01)
#    ax1.set_xlabel('x')
#    ax1.set_ylabel('y')
#    ax2.set_xlabel('x')
#    ax2.set_ylabel('z')
#    ax3.set_xlabel('y')
#    ax3.set_ylabel('z')
#    ax1.set_xlim([-3,3])
#    ax1.set_ylim([-3,3])
#    ax2.set_xlim([-3,3])
#    ax2.set_ylim([0,6])
#    ax3.set_xlim([-3,3])
#    ax3.set_ylim([0,6])
#    ax2.set_title('a={0:.3f}, z={1:.3f}, nH={2:.2f}-{3:.2f}'.format(a,1./a-1,loN,hiN))
#    fig.tight_layout()
#    if label!='':
#        s = 'denseCut_a{0:.3f}_{1:s}.png'.format(a,label)
#    else:
#        s = 'denseCut_a{0:.3f}.png'.format(a)
#    fig.savefig(s,bbox_inches='tight')
#    plt.close(fig)
#
#
## In[ ]:
#
#df.head()
#
#
## In[ ]:
#
#np.sqrt((16.520344**2 + 42.741887**2)/27.098165**2)
#
#
## In[ ]:
#
#loN, hiN = -5, -4
#a,perp,along,speed =[],[],[],[]
#for key, group in groups:
#    ind = (group['loN']==loN) & (group['hiN']==hiN)
#    a.append(key)
#    perp.append(group['vperpStd'][ind].values)
#    along.append(group['valongStd'][ind].values)
#    speed.append(group['speedStd'][ind].values)
#print(perp[0],along[0],speed[0])
#fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
#ax1.plot(a,perp)
#ax2.plot(a,along)
#ax3.plot(a,speed)
#ax1.set_ylabel('v perp std dev')
#ax2.set_ylabel('v along std dev')
#ax3.set_ylabel('speed std dev')
#ax1.set_ylim([0,60])
#ax2.set_ylim([0,60])
#for ax in [ax1,ax2,ax3]:
#    ax.vlines(0.29,0,60,colors='k',linestyles='dashed')
#    ax.vlines(0.32,0,60,colors='k',linestyles='dashed')
#    ax.vlines(0.36,0,60,colors='r',linestyles='dashed')
#    ax.vlines(0.39,0,60,colors='r',linestyles='dashed')
#
#    
#
#
## In[ ]:
#
#df.head()
#
#
## In[ ]:
#
#groups = df.groupby(['loN','hiN'])
#
#
## In[ ]:
#
#df.columns
#
#
## In[ ]:
#
#
#
