
'''
Plots the phase and flows for high z velas
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D


galNums = range(21,30)
expns = ['0.250','0.280','0.340','0.490']

loT, hiT = 10**4, 10**4.5
loN, hiN = 10**-5, 10**-4.5
binrange = [[-10,2],[2,8]]

baseloc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/a{1:s}/'
basename = 'vela2b-{0:d}_GZa{1:s}.h5'

for galNum in galNums:
    
    print galNum
    fig = plt.figure(figsize=(20,20))

    for row, a in enumerate(expns):

        if galNum==24 and a=='0.490':
            a = '0.450'

        loc = baseloc.format(galNum,a)
        fname = basename.format(galNum,a)

        df = pd.read_hdf(loc+fname, 'data')
        index = ( (df['temperature']>loT) & ( df['temperature']<hiT) &
                  (df['density']>loN) & (df['density']<hiN) )
        d = df[index]
        
        # Create the phase diagram
        H, xedges, yedges = np.histogram2d(np.log10(df['density']), np.log10(df['temperature']), bins=50, range=binrange)        
        H = np.rot90(H)
        H = np.flipud(H)
        H = np.ma.masked_where(H==0,H)
        H = np.log10(H)
        ax = fig.add_subplot(4,4,4*row+1)
        mesh = ax.pcolormesh(xedges,yedges,H)
        ax.set_xlim(binrange[0])
        ax.set_ylim(binrange[1])
        ax.set_xlabel('Density')
        ax.set_ylabel('Temperature')
        cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True)
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel('Counts',rotation=270,fontsize=12)

        # Create the scatter plots
        ax = fig.add_subplot(4,4,4*row+2,projection='3d')
        ax.scatter(d['x'],d['y'],d['z'],marker='o',alpha=0.01,color='green')
        ax.view_init(elev=0,azim=0)
        ax.set_xlabel('x')      
        ax.set_ylabel('y')      
        ax.set_zlabel('z')      

        ax = fig.add_subplot(4,4,4*row+3,projection='3d')
        ax.scatter(d['x'],d['y'],d['z'],marker='o',alpha=0.01,color='green')
        ax.view_init(elev=0,azim=90)
        ax.set_xlabel('x')      
        ax.set_ylabel('y')      
        ax.set_zlabel('z')      
        ax.set_title(a)

        ax = fig.add_subplot(4,4,4*row+4,projection='3d')
        ax.scatter(d['x'],d['y'],d['z'],marker='o',alpha=0.01,color='green')
        ax.view_init(elev=90,azim=0)
        ax.set_xlabel('x')      
        ax.set_ylabel('y')      
        ax.set_zlabel('z')      


    fig.tight_layout()
    s = 'vela2b-{0:d}_highz_phase_clouds.png'.format(galNum)
    fig.savefig(s,bbox_inches='tight',dpi=300)







