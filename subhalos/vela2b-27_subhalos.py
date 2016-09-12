

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



loc = '/home/jacob/research/velas/vela2b/vela27/subhalos/'
basename = '{0:s}halos_{1:.3f}.txt'
expns = np.arange(0.200,0.500,0.01)
reds = [1./a - 1 for a in expns]

fig, ax = plt.subplots(1,1,figsize=(5,5))
r = []

xpoint = []
ypoint = []
cpoint = []

outfile = loc+'vela2b-27_cloud1_subhalos.h5'

header = ['a','redshift','MR','x','y','z','r','phi','theta','vr']
#st = '{0:.3f}\t{1:.3f}\t{2:.2f}\t{3:.3f}\t{4:.3f}\t{5:.3f}\t{6:.3f}\t{7:.3f}\t{8:.3f}\n'
st = np.zeros(len(header))
for a,red in zip(expns, reds):
    
    fname = basename.format(loc,a)

    with open(fname,'r') as frot:
        frot.readline()
        frot.readline()
        l = frot.readline().split()
        mhost = float(l[7])
        rhost = float(l[8])

    mvir, rvir, xc, yc, zc, vx, vy, vz = np.loadtxt(fname,skiprows=3,usecols=(7,8,15,16,17,4,5,6), unpack=True)

    r.append(rhost)
    d = [np.sqrt(x**2+y**2+z**2)/rhost for x,y,z in zip(xc,yc,zc)]
    expn = [a for i in d]    
    mass = [m/mhost for m in mvir]

    for i in range(len(mass)):
        if mass[i]>0.1:
            xpoint.append(red)
            ypoint.append(d[i])
            if mass[i]<0.3:
                cpoint.append('black')
            else:
                theta = np.degrees(np.arctan(yc[i]/xc[i]))
                phi = np.degrees(np.arccos(zc[i]/d[i]))
                vr = (vx[i]*xc[i] + vy[i]*yc[i] + vz[i]*zc[i]) / d[i]
                print('Major halo at z = {0:.2f}:\n\tMR = {1:.2f}'.format(red,mass[i]))
                print('\tx = {0:.3f}\ty = {1:.3f}\tz = {2:.3f}'.format(xc[i],yc[i],zc[i]))
                print('\tr = {0:.3f}\tphi = {1:.2f}\ttheta = {2:.2f}'.format(d[i],phi,theta))
                if xc[i]<0 and zc[i]>0:
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
                    st = np.vstack((st,halo))
                    print('\tIn Cloud 1')
                    #f.write(st.format(a,red,mass[i],xc[i],yc[i],zc[i],d[i],phi,theta))
                if mass[i]>0.3 and mass[i]<0.5:
                    cpoint.append('cyan')
                else:
                    cpoint.append('red')

s = ax.scatter(xpoint,ypoint,c=cpoint,marker='o',edgecolor='',s=10)
#ax.plot(expns,r)
ax.set_xlabel('z')
ax.set_ylabel('Distance [Rvir]')


#ax2 = ax.twiny()
#ax2.set_xlim(ax.get_xlim())
#newTickLoc = np.arange(1.0,4.0,0.5)
#aLoc = [1/(z+1.) for z in newTickLoc]
#ax2.set_xticks(aLoc)
#ax2.set_xticklabels(newTickLoc)



#cb = plt.colorbar(s,ax=ax)
#cb.set_label('Halo Mass Ratio')

s = 'vela2b-27_subalos.png'
fig.savefig(s, bbox_inches='tight', dpi=300)

st = np.delete(st, (0), axis=0)
df = pd.DataFrame(st,columns=header)
df.to_hdf(outfile, 'data', mode='w')


