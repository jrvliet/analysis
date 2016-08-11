
from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np

baseloc = '/mnt/cluster/abs/cgm/vela2b/'
baseloc = '/home/jacob/research/velas/vela2b/'
galNum = 27
expn = '0.490'
elev = [0,0,0,90,90,90,180,180,180]
azim = [0,90,180,0,90,180,0,90,180]

fname = '{0:s}vela{1:d}/a{2:s}/vela2b-{1:d}_GZa{2:s}.h5'.format(baseloc,galNum,expn)
df = pd.read_hdf(fname, 'data')

loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4

index = ( (df['temperature']>loT) & (df['temperature']<hiT) &
          (df['density']>loN) & (df['density']<hiN))
d = df[index]
print('Cloud Selected')


fig = plt.figure(figsize=(15,15))
for i in range(0,9):
    print('Plotting with elev={0:d}, azim={1:d}'.format(elev[i],azim[i]))
    
    ax = fig.add_subplot(3,3,i+1,projection='3d')
    ax.scatter(d['x'],d['y'],d['z'],marker='o',alpha=0.01,color='green')
    ax.view_init(elev=elev[i], azim=azim[i])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Elev={0:d}, Azim={1:d}'.format(elev[i],azim[i]))

fig.tight_layout()
fig.savefig('vela2b-{0:d}_viewingAngles.png'.format(galNum),bbox_inches='tight',dpi=300)



