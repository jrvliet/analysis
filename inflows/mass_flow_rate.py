
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


mH = 1.6737e-24
pc2cm = 3.086e18
pc2km = 3.086e13
mSun = 1.989e33


loc = '/home/jacob/research/velas/vela2b/vela27/'
fname = loc+'a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'
rotmat = loc+'a{0:.3f}/rotmat_a{0:.3f}.txt'

times = pd.read_csv('lookback_time.csv')
timeSteps = times['Time until next [Myr]']

expns = np.arange(0.200,0.500,0.01)

aM = []
aR = []

for i,a in enumerate(expns):



    with open(rotmat.format(a), 'r') as f:

        f.readline()
        rvir = float(f.readline().split()[3])


    df = pd.read_hdf(fname.format(a),'data')

    # Calculate mass of each cell in units of Msun
    df['mass'] = df['density']*mH* (df['cell_size']*pc2cm)**3 / mSun 
    
    # Calculate galactocentric distance in units of pc
    df['r'] = np.sqrt( df['x']**2 + df['y']**2 + df['z']**2 )

    # Calculate radial velocity in units of km/s
    df['vr'] = (df['vx']*df['x'] + df['vy']*df['y'] + df['vz']*df['z']) / df['r']

    # Calculate the time between this snapshot and the next in units of s
    timeyr = timeSteps.iloc[i]*1e6
    time = timeyr*3.154e7

    # Calculate the range of each cell in the radial direction in units of pc
    df['range'] = df['vr']*time/pc2km
    
    # Select out all cells that will reach the galaxy
    index = ( df['vr']*time < df['r'] )
    accreting = df[index]

    accretingMass = accreting['mass'].sum()
    accretingRate = accretingMass/timeyr
    
    aM.append(accretingMass)
    aR.append(accretingRate)
    print(a,accretingMass,accretingRate)


print(len(expns),len(aM))
print(type(expns),type(aM))
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))

ax1.plot(expns,aM)
ax1.set_xlabel('a')
ax1.set_ylabel('Accreting Mass [Msun]')
ax1.set_yscale('log')

ax2.plot(expns,aR)
ax2.set_xlabel('a')
ax2.set_ylabel('Accreting Rate [Msun/yr]')

s = 'accretingMass.png'
fig.savefig(s, bbox_inches='tight', dpi=300)




