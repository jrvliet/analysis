
# coding: utf-8

from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
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
    df['vr'] = np.sqrt(df['vxRot']**2 + df['vyRot']**2 + df['vzRot']**2)
    df['vr'] = df['vr']/(kpc2km*rvir)*s2yr
    df['mass'] = df['density']*mH*(df['cell_size']*pc2cm)**3 / mSun
    df['mass'].describe()
    get_ipython().magic(u'timeit')
    shellFlux = np.zeros(len(rPoints))
    print(rPoints)
    for i,dist in enumerate(rPoints):
        shellIndex = (df['rMod']<dist+dr) & (df['rMod']>dist-dr)
        shell = df[shellIndex]
        shellFlux[i] = (shell['mass']*shell['vr'] / dr).sum()
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    ax.plot(rPoints,shellFlux)
    ax.set_xlabel('Distance from Galaxy [Rvir]')
    ax.set_ylabel('Mass Flux [Msun/yr]')
    ax.set_title('a={0:.3f}, z={0:.3f}'.format(a,1./a-1))
    
    print(shellFlux.max())

    s = 'vela2b-27_massFlux_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)



