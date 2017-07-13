
'''
Plots the typical OVI profile
Defined as the mean value in each pixel across 
all z~1 vela snapshots
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


baseloc = '/mnt/cluster/abs/cgm/vela2b/'
subloc = 'vela{0:d}/a0.490/i90/OVI/'

sysabsName = 'vela2b-{0:d}.OVI.los{1:04d}.sysabs'
specName = 'vela2b-{0:d}.OVI.los{1:04d}.OVI{2:d}.spec'

galNums = range(21,30)
losnums = range(1,1000)
trans = [1032,1038]

for tran in trans:
    print('Transition: {0:d}'.format(tran))
    flux = []
    sigma = []
    lengths = []
    i = 0

    for galNum in galNums:

        print('\tGalaxy: {0:d}'.format(galNum))
        loc = baseloc+subloc.format(galNum)

        for losnum in losnums:

            syName = loc+sysabsName.format(galNum,losnum)
            spName = loc+specName.format(galNum,losnum,tran)

            try:
                fi = open(syName,'r')
                fi.close()
                vel,f,s = np.loadtxt(spName,usecols=(1,2,3),unpack=True)
            except IOError:
                pass

            if len(f)<900 and losnum==1:
                print(galNum,losnum,len(f))
            flux.append(f)
            sigma.append(s)
            lengths.append(len(f))
            i += 1

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

    fig.savefig('meanOVI{0:d}profile.png'.format(tran),
                bbox_inches='tight',dpi=300)
    plt.close()
            
            
            


