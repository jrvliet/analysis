

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



loc = '/home/jacob/research/velas/vela2b/vela27/subhalos/'
loc = '/mnt/cluster/abs/Simulations/vela2.1/'
subloc = 'VELA{0:d}/output/ana/'
filename = 'halos_{0:.3f}.txt'

galNums = range(21,30)
expns = np.arange(0.200,0.550,0.01)
reds = [1./a - 1 for a in expns]

fig, ax = plt.subplots(1,1,figsize=(5,5))
r = []

xpoint = []
ypoint = []
cpoint = []


#st = '{0:.3f}\t{1:.3f}\t{2:.2f}\t{3:.3f}\t{4:.3f}\t{5:.3f}\t{6:.3f}\t{7:.3f}\t{8:.3f}\n'

for galNum in galNums:

    outfile = 'vela2b-{0:d}_sigSubHalos.h5'.format(galNum)
    header = ['a','redshift','MR','x','y','z','r','phi','theta','vr','rvir']
    st = np.zeros(len(header))

    for a,red in zip(expns, reds):
        
        fname = loc+subloc.format(galNum)+filename.format(a)

        with open(fname,'r') as frot:
            frot.readline()
            frot.readline()
            l = frot.readline().split()
            mhost = float(l[7])
            rhost = float(l[8])

        mvir, rvir, xc, yc, zc, vx, vy, vz = np.loadtxt(fname,skiprows=3,
                                                usecols=(7,8,15,16,17,4,5,6), 
                                                unpack=True)

        r.append(rhost)
        d = [np.sqrt(x**2+y**2+z**2)/rhost for x,y,z in zip(xc,yc,zc)]
        expn = [a for i in d]    
        mass = [m/mhost for m in mvir]

        for i in range(len(mass)):
            if mass[i]>0.1:
                xpoint.append(red)
                ypoint.append(d[i])
                if mass[i]<0.1:
                    cpoint.append('black')
                else:
                    theta = np.degrees(np.arctan(yc[i]/xc[i]))
                    phi = np.degrees(np.arccos(zc[i]/d[i]))
                    vr = (vx[i]*xc[i] + vy[i]*yc[i] + vz[i]*zc[i]) / d[i]
                    halo = np.zeros(len(header))
                    halo[0] = a
                    halo[1] = z
                    halo[2] = mass[i]
                    halo[3] = xc[i]
                    halo[4] = yc[i]
                    halo[5] = zc[i]
                    halo[6] = d[i]
                    halo[7] = phi
                    halo[8] = theta
                    halo[9] = vr
                    halo[10] = rvir[i]
                    st = np.vstack((st,halo))

    st = np.delete(st, (0), axis=0)
    df = pd.DataFrame(st,columns=header)
    df.to_hdf(outfile, 'data', mode='w')
    print(outfile)


