
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import sqrt
import json
import sys


plot = 0
numSubHalos = 20

loc = '/home/jacob/research/velas/vela2b/vela28/a0.490/i90/'
absfile = 'vela2b-28.0.490.{0:s}.i90.abs_cells.dat'
gasfile = 'vela2b-28_GZa0.490.txt'
ions = ['HI', 'MgII', 'CIV', 'OVI']

outfile = 'vela2b-28.0.490.{0:s}.i90.full_abscells.dat'
outs = '{0:d}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{5:.6f}\t{6:.6f}\t{7:s}\n'
header = 'CellID\tx\t\ty\t\tz\t\tnH\t\tTemp\t\tSNII\n'
print 'Reading in gasfile'
#xLoc, yLoc, zLoc, temp = np.loadtxt(loc+gasfile, skiprows=2, usecols=(1,2,3,8), unpack=True)

with open('cellLocs.json') as f:
    data = json.load(f)

xLoc, yLoc, zLoc = [], [], []
for i in range(len(data[0])):
    xLoc.append(data[0][i])
    yLoc.append(data[1][i])
    zLoc.append(data[2][i])

print len(xLoc)
u = np.linspace(0, 2*np.pi, 100)
v = np.linspace(0, np.pi, 100)


# Open the list of halos
halofile = '{0:s}/../input_0.490.txt'.format(loc)
print 'Reading in halofile'
xhaloraw, yhaloraw, zhaloraw, rvir = np.loadtxt(halofile, usecols=(0,1,2,4), unpack=True)
xhalo, yhalo, zhalo = [], [], []
halonum, rvirhalo = [], []
for i in range(1,numSubHalos):
    xhalo.append( (xhaloraw[i]-xhaloraw[0])*1000. )
    yhalo.append( (yhaloraw[i]-yhaloraw[0])*1000. )
    zhalo.append( (zhaloraw[i]-zhaloraw[0])*1000. )
    halonum.append( i )
    rvirhalo.append( rvir[i] )


for ion in ions:

    
    print '\t'+ion
    fo = open(outfile.format(ion), 'w')
    fo.write(header)
    fname = '{0:s}/{1:s}/{2:s}'.format(loc,ion,absfile.format(ion))

    subhalocounts = 0

#    cellid, nH, t = np.loadtxt(fname, skiprows=1, usecols=(2,7,8), unpack=True)
    print '\tReading in abscell'
    cells = np.loadtxt(fname, skiprows=1)

    t, nH, snII, x, y, z = [], [], [], [], [], []
    inSub = []
    for cell in cells:

        cellid = int(cell[2])
        index = cellid-1
        nH.append(cell[7])
        t.append(cell[8])
        snII.append(cell[9])
        x.append(xLoc[index])
        y.append(yLoc[index])
        z.append(zLoc[index])
        
        # Determine if the cell is part of a subhalo
        subs = []
        inHalo = 0
        for i in range(len(xhalo)):
            dist = sqrt( (xLoc[index]-xhalo[i])**2 + (yLoc[index]-yhalo[i])**2 + (zLoc[index]-zhalo[i])**2 )
            if dist<rvirhalo[i]:
                subs.append(i)
                inHalo = 1
        if inHalo==1:
            subhalocounts += 1
            inSub.append(1)
        else:
            inSub.append(0)
        
            
        if len(subs)>0:
            s = ''
            for i in range(len(subs)):
                s += '{0:d}\t'.format(subs[i])
        else:
            s = '0'
        fo.write(outs.format(cellid,xLoc[index],yLoc[index],zLoc[index],cell[7],cell[8],cell[9],s))

    maxpoint = max(x)
    minpoint = min(x)

    if max(y)>maxpoint:
        maxpoint = max(y)
    if max(z)>maxpoint:
        maxpoint = max(z)

    if min(y)<minpoint:
        minpoint = min(y)
    if min(z)<minpoint:
        minpoint = min(z)

    if plot==1:
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(x, y, z, c=inSub, s=20, marker='o', cmap='gray')

        minpoint = -300
        maxpoint = 300
        ax.set_xlim([minpoint,maxpoint])
        ax.set_ylim([minpoint,maxpoint])
        ax.set_zlim([minpoint,maxpoint])

        # Plot subhalos
        for i in range(0,len(rvirhalo)):
            xs = rvirhalo[i] * np.outer(np.cos(u), np.sin(v)) + xhalo[i]
            ys = rvirhalo[i] * np.outer(np.sin(u), np.sin(v)) + yhalo[i]
            zs = rvirhalo[i] * np.outer(np.ones(np.size(u)), np.cos(v)) + zhalo[i]

            ax.plot_surface(xs,ys,zs, rstride=1, cstride=1, alpha=0.9, color='r', linewidth=0)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_aspect('equal')

        plt.show()
        sys.exit()
        

    print '\t\t{0:d} of {1:d} ({2:.2%}) cells in subhalos'.format(subhalocounts, len(cells), subhalocounts/float(len(cells)))
    fo.close()

