# coding: utf-8 
from __future__ import print_function 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
pd.options.mode.chained_assignment = None 

kpc2km = 3.086e16 
pc2cm = 3.086e18 
mH = 1.6737e-24 
mSun = 1.989e33 
s2yr = 3.154e7 

dataloc = '/home/jacob/research/velas/vela2b/vela27/' 
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5' 
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expns = np.arange(0.200,0.500,0.01)
dr = 0.01
rPoints = np.arange(0,5,0.1)

eras = 'Formation Merger Starburst Relax'.split()
eraSplits = [0.29, 0.33, 0.36]
eraBinning = np.digitize(expns,eraSplits)

tempLabels = 'cool warm hot'.split()
tempEdges = np.linspace(3.5,4.5,4)


for i,a in enumerate(expns):
    fname = dataloc + filename.format(a)
    rname = dataloc + rotname.format(a)

    df = pd.read_hdf(fname,'data')

    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    df['rMod'] = df['r']/rvir
    #df['rkm'] = df['r']*kpc2km
    df['vr'] = (df['vxRot']*df['xRot'] + df['vyRot']*df['yRot'] +
                        df['vzRot']*df['zRot']) / df['r']
    df['vr'] = df['vr']/(kpc2km*rvir)*s2yr
    df['mass'] = df['density']*mH*(df['cell_size']*pc2cm)**3 / mSun
    df['mass'].describe()

    
    # Select out the filament space regeme
    spaceIndex = (df['theta']<80)
    
    filamentIn = np.zeros(len(rPoints),len(tempLabels))
    filamentOut = np.zeros(len(rPoints),len(tempLabels))
    shellIn = np.zeros(len(rPoints))
    shellOut = np.zeros(len(rPoints))

    for j,tempLabel in enumerate(tempLabels):
        loT = tempEdges[j]
        hiT = tempEdges[j+1]
        tempIndex = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
 
        for k,dist in enumerate(rPoints):
            shellIndex = (df['rMod']<dist+dr) & (df['rMod']>dist-dr)
            inIndex = (df['vr']<0)
            outIndex = (df['vr']>0)

            shellInfall = df[shellIndex & inIndex]
            shellOutflow = df[shellIndex & outIndex]

            filamentInfall = df[shellIndex & inIndex & tempIndex & spaceIndex]
            filamentOutflow = df[shellIndex & outIndex & tempIndex & spaceIndex]
            
            shellIn[k] = (shellInfall['mass']*shellInfall['vr'] / dr).sum()
            shellOut[k] = (shellOutflow['mass']*shellOutflow['vr'] / dr).sum()
    
            filamentIn[k] = (filamentInfall['mass']*filamentInfall['vr'] / dr).sum()
            filamentOut[k] = (filamentOutflow['mass']*filamentOutflow['vr'] / dr).sum()
    
    #print('In:  {0:.3f} - {1:.3f}'.format(shellIn.min(),shellIn.max()), end='\t')
    #print('Out: {0:.3f} - {1:.3f}'.format(shellOut.min(),shellOut.max()))

    fig,(ax1,ax2,ax3,ax4) = plt.subplots(1,4,figsize=(20,5))
    ax1.plot(rPoints,np.abs(shellIn),color='blue',label='Inflow')
    ax1.plot(rPoints,shellOut,color='red',label='Outflow')
    ax1.set_title('All')

    for j, (ax,tempLabel) in enumerate(zip([ax2,ax3,ax4],tempLabel)):
        ax1.plot(rPoints,np.abs(filamentIn[:,j]),color='blue',label='Inflow')
        ax1.plot(rPoints,filamentOut[:,j],color='red',label='Outflow')
        ax.set_title('{0:s} Filament'.format(tempLabel.title()))
    
    #for ax in [ax1,ax2]:
    for ax in [ax]:
        ax.set_xlabel('Distance from Galaxy [Rvir]')
        ax.set_ylabel('Mass Flux [Msun/yr]')
        ax.set_xlim([0,5])
        ax.set_ylim([0,200])
    ax4.legend(loc='upper right')

    #ax1.set_title('Inflow')
    #ax2.set_title('Outflow')
    fig.suptitle('a={0:.3f}, z={1:.3f}'.format(a,1./a-1))
    fig.tight_layout()
    
    #ax.set_title('a={0:.3f}, z={1:.3f}'.format(a,1./a-1))

    ax.text(0.75,0.5,eras[eraBinning[i]],transform=ax.transAxes)
    s = 'vela2b-27_massFlux_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


