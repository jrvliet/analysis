

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
for i,(field,label) in enumerate(zip(fields,labels)):
    
    print(field)
    print(type(cells))
    print(type(cells[field]))
    print(cells[field].head())

    # First row is entire sample
    ax[0,i].hist(cells[field],bins=30,normed=True,histtype='step',color='black',label='Full')

    # Second row is major/minor split as defined by cells
    ax[1,i].hist(cellMajor[field],bins=30,normed=True,histtype='step',
            color='red',label='Cell Major')
    ax[1,i].hist(cellMinor[field],bins=30,normed=True,histtype='step',
            color='blue',label='Cell Minor')

    # Third row is major/minor split as defined by LOS
    ax[2,i].hist(losMajor[field],bins=30,normed=True,histtype='step',
            color='red',label='LOS Major')
    ax[2,i].hist(losMinor[field],bins=30,normed=True,histtype='step',
            color='blue',label='LOS Minor')

# Label
for a in ax.flatten():
    a.set_xlabel('{0:s} Velocity [km/s]'.format(label))
    a.set_xlim([-800,800])
    a.legend(loc='best')

fig.tight_layout()
fig.savefig('oviSphericalVelocitiesDist.png',bbox_inches='tight',dpi=300)
        
    





