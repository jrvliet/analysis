
'''
Plots the typical OVI profile
Defined as the mean value in each pixel across 
all z~1 vela snapshots
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def azimuthal(phi):
    if phi<90:
        return phi
    elif phi<180:
        return 180.-phi
    elif phi<270:
        return phi-180.
    else:
        return 360.-phi

baseloc = '/mnt/cluster/abs/cgm/vela2b/'
subloc = 'vela{0:d}/a0.490/i90/OVI/'

sysabsName = 'vela2b-{0:d}.OVI.los{1:04d}.sysabs'
specName = 'vela2b-{0:d}.OVI.los{1:04d}.OVI{2:d}.spec'
linesHeader = 'losnum D phi incline'.split()

galNums = range(21,30)
losnums = range(1,1000)
trans = [1032,1038]

azCutoff = 30.
loD,hiD = 20,250

for tran in trans:
    print('Transition: {0:d}'.format(tran))
    majorFlux,minorFlux = [],[]
    majorSigma,minorSigma = [],[]
    i = 0

    for galNum in galNums:

        print('\tGalaxy: {0:d}'.format(galNum))
        loc = baseloc+subloc.format(galNum)
        lines = loc+'lines.info'
        try:
            lines = pd.read_csv(lines,sep='\s+',skiprows=2,names=linesHeader)
            lines['az'] = lines['phi'].apply(azimuthal)
            print(galNum)
        except IOError:
            continue

        for losnum in losnums:

            impact = lines['D'].iloc[losnum-1]
            if impact>=loD and impact<=hiD:
                syName = loc+sysabsName.format(galNum,losnum)
                spName = loc+specName.format(galNum,losnum,tran)

                try:
                    fi = open(syName,'r')
                    fi.close()
                    vel,f,s = np.loadtxt(spName,usecols=(1,2,3),unpack=True)
                
                    az = lines['az'].iloc[losnum-1]
                    if az<azCutoff:
                        majorFlux.append(f)
                        majorSigma.append(s)
                    else:
                        minorFlux.append(f)
                        minorSigma.append(s)
                        
                except IOError:
                    continue

            i += 1

    
    for flux,sigma,name in zip([majorFlux,minorFlux],
                               [majorSigma,minorSigma],
                               ['Major','Minor']):

        flux = np.asarray(flux)
        sigma = np.asarray(sigma)
        top = np.sum(flux/sigma**2,axis=0)
        bot = np.sum(sigma**(-2),axis=0)
        y = top/bot
        yerr = 1./np.sqrt(bot**2)
        #y = np.mean(flux,axis=0)
        #yerr = np.std(flux,axis=0)/np.sqrt(tot.shape[0])

        fig,ax = plt.subplots(1,1,figsize=(5,5))
        ax.errorbar(vel,y,yerr=yerr)
        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel('Flux')
        ax.set_ylim([0,1.2])

        fig.savefig('meanOVI{0:d}_{1:s}_Limitsprofile.png'.format(tran,name),
                    bbox_inches='tight',dpi=300)
        plt.close()
            
            
            


