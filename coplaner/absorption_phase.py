'''
Plots the phase of gas based on the geometry of the situation
pulled from the gas box

coplaner = 20 kpc thick cylinder centered on galaxy, inf radius
outflows = 4o kpc radius cylinder cenetered on galaxy, inf thickness

'''


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

def select_regions(x, y, z, aOutflow, aPlane, rmin, rmax): 

    # To be in the plane, cell must not be in the cone with 
    # theta = arctan(aPlane)
    planeInds = ( (x**2+y**2)/aPlane > z**2 &
                  (x**2 + y**2 + z**2) > rmin &
                  (x**2 + y**2 + z**2) < rmax )

    # To be in the in outflows, cell must be in the cone with
    # theta = arctan(aOutflow)
    outflowInds = ( (x**2+y**2)/aOutflow < z**2 &
                    (x**2 + y**2 + z**2) > rmin &
                    (x**2 + y**2 + z**2) < rmax )


    # To be in the void, cell must be in the plane cone and
    # not in the outflow cone
    voidInds = ( ~planeInds & ~outflowInds &
                    (x**2 + y**2 + z**2) > rmin &
                    (x**2 + y**2 + z**2) < rmax )

    return planeInds, outflowInds, voidInds

def plot_hist( dense, temp, ind, ax, title ):
    
    numbins = 50
    binrange = [[-10,3],[2,8]]

    n = dense[ind]
    t = temp[ind]

    h, xedges, yedges = np.histogram2d(n, t, bins=numbins, range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
    mesh = ax.pcolormesh(xedges, yedges, h)
    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    ax.set_title(title)



coplaneZCut = 20    # [kpc]
outflowRCut = 40    # [kpc]
galRCut = 10        # [kpc]
galZCut = 5         # [kpc] 

# Define the cones
outflowAngle = 45
planeAngle = 70
aOut = np.tan(np.radians(outflowAngle))
aPlane = np.tan(np.radians(planeAngle))

# Read in the gas box
print 'Read in gas box'
galID = 28
expn = '0.490'
gasbox = 'vela2b-{0:d}_GZa{1:s}.h5'.format(galID, expn)
d = pd.read_hdf(gasbox, 'data')

ions = ['MgII', 'CIV', 'OVI']
absfile = './i{0:d}/{1:s}/vela2b-{2:d}.{3:s}.{1:s}.abs_cells.h5'

xbox = d['x']
ybox = d['y']
zbox = d['z']
dense = np.log10(d['density'])
temp = np.log10(d['temperature'])

# Read in the rotation matrix
print 'Read in rotation matrix'
rotmat = np.zeros((3,3))
rotmatFile = 'rotmat_a{0:s}.txt'.format(expn)
with open(rotmatFile, 'r') as f:
    f.readline()
    l = f.readline().split()
    rvir = float(l[3])
    a11 = float(l[4])
    a12 = float(l[5])
    a13 = float(l[6])
    a21 = float(l[7])
    a22 = float(l[8])
    a23 = float(l[9])
    a31 = float(l[10])
    a32 = float(l[11])
    a33 = float(l[12])


# This matrix describes the rotation from the box frame
# to the galaxy frame
# Convert all the cell locations from the box frame to the 
# galaxy frame
readin = 1
print 'Convert to galaxy frame'
if readin == 0:
    xgal, ygal, zgal = [], [], []
    f = open('galaxy.coords', 'w')
    s = '{0:f}\t{1:f}\t{2:f}\n'
    for i in range(len(xbox)):
        x = a11*xbox[i] + a12*ybox[i] + a13*zbox[i] 
        y = a21*xbox[i] + a22*ybox[i] + a23*zbox[i] 
        z = a31*xbox[i] + a32*ybox[i] + a33*zbox[i] 
        f.write(s.format(x, y, z))
        
        xgal.append(x)
        ygal.append(y)
        zgal.append(z)
    f.close()
    x = np.array(xgal)
    y = np.array(ygal)
    z = np.array(zgal)
else:
    x, y, z = np.loadtxt('galaxy.coords', usecols=(0,1,2), unpack=True)




rmin = 0.2*rvir
rmax = 0.75*rvir

fig, axes = plt.subplots(4,4, figsize=(16,16))


axs = axes[0]    

    # Plot all gas in the box
    planeInds, outflowInds, voidInds = select_regions(x, y, z, aOutflow, aPlane, rmin, rmax) 
    allInds = np.array( [True for i in range(len(x))] )
    plot_hist( dense, temp, planeInds, axs[0], 'All Coplanar Gas') 
    plot_hist( dense, temp, outflowInds, axs[1], 'All Outflow Gas') 
    plot_hist( dense, temp, voidInds, axs[2], 'All Void Gas') 
    plot_hist( dense, temp, allInds, axs[3], 'All Gas') 
    
    
plot_hist( dense, temp, ind, ax, title ):
# Plot coplanar gas
hPlane, xedges, yedges = np.histogram2d(nPlane, tPlane, bins=numbins, range=binrange)
h = np.rot90(hPlane)
h = np.flipud(h)
h = np.ma.masked_where(h==0,h)
h = np.log10(h)
mesh = ax1.pcolormesh(xedges, yedges, h)
cbar = plt.colorbar(mesh, ax=ax1, use_gridspec=True)
ax1.set_title('Coplanar Gas')

# Plot gas in outflow region
hOutflow, xedges, yedges = np.histogram2d(nOutflow, tOutflow, bins=numbins, range=binrange)
h = np.rot90(hOutflow)
h = np.flipud(h)
h = np.ma.masked_where(h==0,h)
h = np.log10(h)
mesh = ax2.pcolormesh(xedges, yedges, h)
cbar = plt.colorbar(mesh, ax=ax2, use_gridspec=True)
ax2.set_title('Outflow Gas')

# Plot gas not in the two regions
hVoid, xedges, yedges = np.histogram2d(nVoid, tVoid, bins=numbins, range=binrange)
h = np.rot90(hVoid)
h = np.flipud(h)
h = np.ma.masked_where(h==0,h)
h = np.log10(h)
mesh = ax3.pcolormesh(xedges, yedges, h)
cbar = plt.colorbar(mesh, ax=ax3, use_gridspec=True)
ax3.set_title('Void Gas')

# Plot all gas
hAll, xedges, yedges = np.histogram2d(dense, temp, bins=numbins, range=binrange)
h = np.rot90(hAll)
h = np.flipud(h)
h = np.ma.masked_where(h==0,h)
h = np.log10(h)
mesh = ax4.pcolormesh(xedges, yedges, h)
cbar = plt.colorbar(mesh, ax=ax4, use_gridspec=True)
ax4.set_title('All Gas')

# Plot Normalized coplane
h = hPlane / hAll
h[np.isnan(h)] = 0.0
h = np.rot90(h)
h = np.flipud(h)
h = np.ma.masked_where(h==0,h)
mesh = ax5.pcolormesh(xedges, yedges, h)
cbar = plt.colorbar(mesh, ax=ax5, use_gridspec=True)
ax5.set_title('Normed Plane Gas')


# Plot normalized outflows
h = hOutflow / hAll
h[np.isnan(h)] = 0.0
h = np.rot90(h)
h = np.flipud(h)
h = np.ma.masked_where(h==0,h)
mesh = ax6.pcolormesh(xedges, yedges, h)
cbar = plt.colorbar(mesh, ax=ax6, use_gridspec=True)
ax6.set_title('Normed Outflow Gas')


# Plot normalized void
h = hVoid / hAll
h[np.isnan(h)] = 0.0
h = np.rot90(h)
h = np.flipud(h)
h = np.ma.masked_where(h==0,h)
mesh = ax7.pcolormesh(xedges, yedges, h)
cbar = plt.colorbar(mesh, ax=ax7, use_gridspec=True)
ax7.set_title('Normed Void Gas')

# Plot gas that is within both plane and outflow
hBoth, xedges, yedges = np.histogram2d(nBoth, tBoth, bins=numbins, range=binrange)
h = np.rot90(hBoth)
h = np.flipud(h)
h = np.ma.masked_where(h==0,h)
h = np.log10(h)
mesh = ax8.pcolormesh(xedges, yedges, h)
cbar = plt.colorbar(mesh, ax=ax8, use_gridspec=True)
ax8.set_title('Outflow and Plane Gas')


for axs in axes:

    for ax in axs:
        ax.set_xlim(binrange[0])
        ax.set_ylim(binrange[1])
        ax.set_xlabel('log nH')
        ax.set_ylabel('log T')



plt.tight_layout()
s = 'inflow_outflow_phase.png'
plt.savefig(s, dpi=300, bbox_inches='tight')


