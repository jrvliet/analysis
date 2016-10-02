
'''
Calculates virial velocity and temperature
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


mSun2g = 1.98855e33
kpc2cm = 3.086e21
G = 6.6743e-8
mp = 1.6726e-24
kb = 1.3807e-16
cms2kms = 1e-5


loc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expns = np.arange(0.200,0.500,0.01)

outfile = 'vela2b-27_virialProps.dat'
header = ['a','z','mvir','rvir','vvir','tvir']
vir = np.zeros(len(header))

for a in expns:
    

    fname = loc+filename.format(a)
    with open(fname,'r') as f:
        f.readline()
        l = f.readline().split()
        z = float(l[1])
        mvir = float(l[2])
        rvir = float(l[3])

    mvir_g = mvir*mSun2g
    rvir_cm = rvir*kpc2cm
    vVir_cms = np.sqrt( G*mvir_g/rvir_cm )
    tVir = mp*vVir_cms**2/(2*kb)
    
    vVir = vVir_cms*cms2kms

    halo = np.zeros(len(header))
    halo[0] = a
    halo[1] = z
    halo[2] = mvir
    halo[3] = rvir
    halo[4] = vVir
    halo[5] = tVir

    vir = np.vstack((vir,halo))

    print('a = {0:.3f}\tMvir = {1:.3e}\tRvir = {2:.3f}\tVvir = {3:.3f}\tTvir = {4:.3e}'.format(a,mvir,rvir,vVir,tVir))
    

vir = np.delete(vir,(0),axis=0)
df = pd.DataFrame(vir,columns=header)
df.to_csv(outfile)

fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))
ax1.plot(df['z'],df['mvir'])
ax2.plot(df['z'],df['rvir'])
ax3.plot(df['z'],df['vvir'])
ax4.plot(df['z'],df['tvir'])

ax1.invert_xaxis()
ax2.invert_xaxis()
ax3.invert_xaxis()
ax4.invert_xaxis()

ax1.set_xlabel('z')
ax2.set_xlabel('z')
ax3.set_xlabel('z')
ax4.set_xlabel('z')

ax1.set_ylabel('Mvir')
ax2.set_ylabel('Rvir')
ax3.set_ylabel('Vvir')
ax4.set_ylabel('Tvir')

ax1.set_yscale('log')
ax4.set_yscale('log')

fig.tight_layout()
fig.savefig('virial_props.png',bbox_inches='tight',dpi=300)









