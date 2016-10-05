
'''
Determines the mass flux through several regions in 
the filament
'''

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


pd.options.mode.chained_assignment = None

# Conversion factors
Myr2s = 3.154e13   # Number of s in Myr
pc2km = 3.086e13   # Number of km in a pc

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZ{0:.3f}.rot.h5'
rotmatfile = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expns0 = range(20,50)
expns = [i/100. for i in expns0]

# Define phase limits
loT, hiT = 10**3.5, 10**4.5
limits = pd.read_csv('cloudLimits.csv')

# Read in the table of time between snapshots
times = pd.read_csv('lookback_time.csv')

numbins = 4

for a in expns:

    # Read in the box
    fname = baseloc+filename.format(a)
    df = pd.read_hdf(fname, 'data')

    # Get rvir
    with open(baseloc+rotmatfile.format(a)) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    # Get density limits
    nRange = limits[limits.expn==a]
    loN = 10**(nrange['loN'].values[0])
    hiN = 10**(nrange['hiN'].values[0])

    # Select out the filament
    tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
    densInd = (df['density']>loN) & (df['density']<hiN)
    spacInd = (df['zRot']>0) & (df['theta']<80)

    cl = df[tempInd & densInd & spacInd]
    
    # Calculate the projected new location of the gas
    # based on current velocity
    timestepMyr = times[times.a==a]['Time until next [Myr]'].values[0]
    timestepS = timestepMyr*Myr2s
    cl['xNewMod'] = ((cl['x']*pc2km + cl['vx']*timestepS)/pc2km)/rvir
    cl['yNewMod'] = ((cl['y']*pc2km + cl['vy']*timestepS)/pc2km)/rvir
    cl['zNewMod'] = ((cl['z']*pc2km + cl['vz']*timestepS)/pc2km)/rvir

    cl['rNewMod'] = np.sqrt(cl['xNewMod']**2 + cl['yNewMod']**2 + cl['zNewMod']**2 )


    # Create the bins
    cl['rmod'] = cl['r']/rvir

    rbinEdges = np.linspace(0,cl['rmod'].max(),numbins)
    
    # Bin by r
    groups = cl.groupby( pd.cut(cl['rmod'],rbinEdges) )

    for name, group in groups:

        
    
    
    
    
    





