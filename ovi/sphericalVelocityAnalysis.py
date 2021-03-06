

'''
Plots distribution of OVI absorbing cells' 
spherical velocities
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

loc = '/home/jacob/research/code/analysis/ovi/'
filename = 'vela2b-{0:d}.a0.490.i90_OVIcellsGalFrame.h5'
galNums = range(21,30)

for i,galNum in enumerate(galNums):

    fname = loc+filename.format(galNum)
    try:
        df = pd.read_hdf(fname,'data')
    except IOError:
        continue
    
    if i==0:
        cells = df.copy()
    else:
        cells = pd.concat([cells,df])

print(len(cells))

cellMajor = cells[cells['cellPlane']]
cellMinor = cells[cells['cellOutflow']]
losMajor = cells[cells['losPlane']]
losMinor = cells[cells['losOutflow']]

fig,ax = plt.subplots(3,3,figsize=(15,15),sharey='row',sharex=True)

fields = 'vr vphi vtheta'.split()
labels = 'Radial Polar Rotating'.split()
bins = np.linspace(-800,800,50)

cols = ('vr_mean vr_std vr_skew'.split() + 
        'vPolar_mean vPolar_std vPolar_skew'.split() + 
        'vRotate_mean vRotate_std vRotate_skew'.split())
prefix = 'vr vPolar vRotate'.split()
rows = 'full cell_major cell_minor los_major los_minor'.split()

stats = pd.DataFrame(columns=cols,index=rows)

for i,(field,label) in enumerate(zip(fields,labels)):
    
    # First row is entire sample
    ax[0,i].hist(cells[field],bins=bins,normed=True,histtype='step',color='black',label='Full')

    stats[prefix[i]+'_mean'].loc['full'] = cells[field].mean()
    stats[prefix[i]+'_std'].loc['full'] = cells[field].std()
    stats[prefix[i]+'_skew'].loc['full'] = cells[field].skew()
    print('Field: {0:s}\t Full\t Num Nan: {1:d}'.format(field, 
            cells[field].isnull().sum()))

    # Second row is major/minor split as defined by cells
    ax[1,i].hist(cellMajor[field],bins=bins,normed=True,histtype='step',
            color='red',label='Cell Major')
    ax[1,i].hist(cellMinor[field],bins=bins,normed=True,histtype='step',
            color='blue',label='Cell Minor')

    stats[prefix[i]+'_mean'].loc['cell_major'] = cellMajor[field].mean()
    stats[prefix[i]+'_std'].loc['cell_major'] = cellMajor[field].std()
    stats[prefix[i]+'_skew'].loc['cell_major'] = cellMajor[field].skew()

    stats[prefix[i]+'_mean'].loc['cell_minor'] = cellMinor[field].mean()
    stats[prefix[i]+'_std'].loc['cell_minor'] = cellMinor[field].std()
    stats[prefix[i]+'_skew'].loc['cell_minor'] = cellMinor[field].skew()

    print('Field: {0:s}\t Cell Major\t Num Nan: {1:d}'.format(field, 
            cellMajor[field].isnull().sum()))
    print('Field: {0:s}\t Cell Minor\t Num Nan: {1:d}'.format(field, 
            cellMinor[field].isnull().sum()))

    # Third row is major/minor split as defined by LOS
    ax[2,i].hist(losMajor[field],bins=bins,normed=True,histtype='step',
            color='red',label='LOS Major')
    ax[2,i].hist(losMinor[field],bins=bins,normed=True,histtype='step',
            color='blue',label='LOS Minor')

    stats[prefix[i]+'_mean'].loc['los_major'] = losMajor[field].mean()
    stats[prefix[i]+'_std'].loc['los_major'] = losMajor[field].std()
    stats[prefix[i]+'_skew'].loc['los_major'] = losMajor[field].skew()

    stats[prefix[i]+'_mean'].loc['los_minor'] = losMinor[field].mean()
    stats[prefix[i]+'_std'].loc['los_minor'] = losMinor[field].std()
    stats[prefix[i]+'_skew'].loc['los_minor'] = losMinor[field].skew()

    print('Field: {0:s}\t LOS Major\t Num Nan: {1:d}'.format(field, 
            losMajor[field].isnull().sum()))
    print('Field: {0:s}\t LOS Minor\t Num Nan: {1:d}'.format(field, 
            losMinor[field].isnull().sum()))

# Label
for a in ax.flatten():
    a.set_xlim([-800,800])
    a.legend(loc='best')

for a,label in zip(ax[-1,:],labels):
    a.set_xlabel('{0:s} Velocity [km/s]'.format(label))

for a in ax[:,0]:
    a.set_ylabel('Relative Counts')

fig.subplots_adjust(hspace=0.05,wspace=0)
fig.savefig('oviSphericalVelocitiesDist.png',bbox_inches='tight',dpi=300)

stats.to_csv('oviSphericalVelocitiesBulk.csv')




