
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as mt


pd.options.mode.chained_assignment = None

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'


expns = np.arange(0.200,0.500,0.01)

minRho,maxRho = 0,3
numRhoBins = 50
rhoBins = np.linspace(minRho,maxRho,numRhoBins+1)

loT,hiT = 3.5,4.5
lowestN,highestN = -5.5,-2.5
numNBins = 6
nBins = np.linspace(lowestN,highestN,numNBins+1)

for a in expns:

    print(a)

    rname = dataloc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname,'data')

    tempInds = (df['temperature']<10**hiT) & (df['temperature']>10**loT)
    spacInds = (df['theta']<80)


    fig,axes = plt.subplots(1,2,figsize=(10,5))
    for j in range(numNBins):
        loN = nBins[j]
        hiN = nBins[j+1]
        densInds = (df['density']<10**hiN) & (df['density']>10**loN)
        fil = df[tempInds & densInds & spacInds]
        fil['rho'] = np.sqrt(fil['xRot']**2 + fil['yRot']**2) / rvir

        r = np.zeros(numRhoBins)
        t = np.zeros(numRhoBins)
        for i in range(numRhoBins):
            minR = rhoBins[i]
            maxR = rhoBins[i+1]
            shell = fil[(fil['rho']>minR) & (fil['rho']<maxR)]
            r[i] = shell['rho'].mean()
            t[i] = np.log10(mt.gmean(shell['temperature']))

        axes[0].plot(r,t,label='nBin = {0:d}'.format(j+1))
        axes[0].set_xlabel('Distance from filament [Rvir]')
        axes[0].set_ylabel('Mean Temperature')
        axes[0].set_xlim([0,3])
        axes[0].set_ylim([loT,hiT])
            
        axes[1].plot(r*rvir,t,label='nBin = {0:d}'.format(j+1))
        axes[1].set_xlabel('Distance from filament [kpc]')
        axes[1].set_ylabel('Mean Temperature')
        #axes[1].set_xlim([0,3])
        axes[1].set_ylim([loT,hiT])
        
    axes[0].legend(loc='best')
    fig.tight_layout()
    s = 'tempProfile_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


