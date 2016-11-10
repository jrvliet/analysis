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

for a in expns:
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

    shellIn = np.zeros(len(rPoints))
    shellOut = np.zeros(len(rPoints))
    
    for i,dist in enumerate(rPoints):
        shellIndex = (df['rMod']<dist+dr) & (df['rMod']>dist-dr)
        inIndex = (df['vr']<0)
        outIndex = (df['vr']>0)

#        print(len(inIndex))
#        print(len(outIndex))
        shellInfall = df[shellIndex & inIndex]
        shellOutflow = df[shellIndex & outIndex]
        
        shellIn[i] = (shellInfall['mass']*shellInfall['vr'] / dr).sum()
        shellOut[i] = (shellOutflow['mass']*shellOutflow['vr'] / dr).sum()
    
    print('In:  {0:.3f} - {1:.3f}'.format(shellIn.min(),shellIn.max()), end='\t')
    print('Out: {0:.3f} - {1:.3f}'.format(shellOut.min(),shellOut.max()))

    #fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
    #ax1.plot(rPoints,np.abs(shellIn))
    #ax2.plot(rPoints,shellOut)

    fig,ax = plt.subplots(1,1,figsize=(5,5))
    ax.plot(rPoints,np.abs(shellIn),color='blue',label='Inflow')
    ax.plot(rPoints,shellOut,color='red',label='Outflow')
    
    #for ax in [ax1,ax2]:
    for ax in [ax]:
        ax.set_xlabel('Distance from Galaxy [Rvir]')
        ax.set_ylabel('Mass Flux [Msun/yr]')
        ax.set_xlim([0,5])
        ax.set_ylim([0,200])

    #ax1.set_title('Inflow')
    #ax2.set_title('Outflow')
    #fig.suptitle('a={0:.3f}, z={1:.3f}'.format(a,1./a-1))
    #fig.tight_layout()
    
    ax.set_title('a={0:.3f}, z={1:.3f}'.format(a,1./a-1))
    ax.legend(loc='upper right')

    s = 'vela2b-27_massFlux_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


