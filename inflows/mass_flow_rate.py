
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


mH = 1.6737e-24
pc2cm = 3.086e18
mSun = 1.989e33


loc = '/home/jacob/research/velas/vela2b/vela27/'
fname = loc+'a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'
rotmat = loc+'a{0:.3f}/rotmat_a{0:.3f}.txt'

times = pd.read_csv('lookback_time.csv')
timeSteps = times['Time until next [Myr]']

expns = np.arange(0.200,0.500,0.01)

aM = []

for i,a in enumerate(expns):

    time = timeSteps.iloc[i]*1e6*3.154e7


    with open(rotmat.format(a), 'r') as f:

        f.readline()
        rvir = float(f.readline().split()[3])


    df = pd.read_hdf(fname.format(a),'data')

    df['mass'] = df['density']*mH* (df['cell_size']*pc2cm)**3 / mSun 
    
    df['r'] = np.sqrt( df['x']**2 + df['y']**2 + df['z']**2 )

    df['vr'] = (df['vx']*df['x'] + df['vy']*df['y'] + df['vz']*df['z']) / df['r']

    df['range'] = df['vr']*time
    
    index = ( df['vr']*time < df['r'] )

    accreting = df[index]

    accretingMass = accreting['mass'].sum()
    
    aM.append(accretingMass)
    print(a,accretingMass)


print(len(expns),len(aM))
print(type(expns),type(aM))
fig, ax = plt.subplots(1,1,figsize=(5,5))

ax.plot(expns,aM)
ax.set_xlabel('a')
ax.set_ylabel('Accreting Mass')
ax.set_yscale('log')

s = 'accretingMass.png'
fig.savefig(s, bbox_inches='tight', dpi=300)




