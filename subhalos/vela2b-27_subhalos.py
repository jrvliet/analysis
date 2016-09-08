

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt



loc = '/home/jacob/research/velas/vela2b/vela27/subhalos/'
basename = '{0:s}halos_{1:.3f}.txt'
expns = np.arange(0.200,0.500,0.01)
red = [1./a - 1 for a in expns]

fig, ax = plt.subplots(1,1,figsize=(5,5))
r = []

for a,z in zip(expns, red):
    
    print('a = {0:.3f}'.format(a))

    fname = basename.format(loc,a)

    with open(fname,'r') as f:
        f.readline()
        f.readline()
        l = f.readline().split()
        mhost = float(l[7])
        rhost = float(l[8])

    mvir, rvir, xc, yc, zc = np.loadtxt(fname,skiprows=3,usecols=(7,8,15,16,17), unpack=True)

    r.append(rhost)
    d = [np.sqrt(x**2+y**2+z**2) for x,y,z in zip(xc,yc,zc)]
    expn = [a for i in d]    
    mass = [m/mhost for m in mvir]

    ax.scatter(expn[:10],d[:10],c=mass[:10], marker='o')

ax.plot(expns,r)
ax.set_xlabel('a')
ax.set_ylabel('Distance [kpc]')

s = 'vela2b-27_subalos.png'
fig.savefig(s, bbox_inches='tight', dpi=300)



