
'''
A code to quantify the clumping of absorbing cells in spatial location, 
density, temperature, and metallicity
'''
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from scipy.stats import kde
from jenks import jenks

def get_delta_s(xen, yen, zen, xex, yex, zex, x, y, z, imp):
    
    '''
    Returns the distance the point given by (x,y,z) is from
    (xen, yen, zen) along the LOS defined by the starting
    and ending points
    '''
    
    # Distnace from galaxy to cell
    r = sqrt( x*x + y*y + z*z )

    # Distance from galaxy to exit point
    end = sqrt( xex**2 + yex**2 + zex**2 )

    # Length of LOS
    l = sqrt( (xen-xex)**2 + (yen-yex)**2 + (zen-zex)**2 )

    # Distance from entry point to cell
#    s = l - (sqrt( r**2 - imp**2 ) + sqrt( end**2 - imp**2) )
    s = sqrt( (xen-x)**2 + (yen-y)**2 + (zen-z)**2 )

    return s/l

testLoc = '/home/hyades/jrvander/exampleData/'

ions = ['HI', 'MgII', 'CIV', 'OVI']

galID = 'vela2b-28'
expn = '0.510'
inc = '10'

fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ((ax31, ax32),(ax33,ax34)) = plt.subplots(2,2)
ax3 = (ax31, ax32, ax33, ax34)
numbins = 20

for ionnum, ion in enumerate(ions):

    print ion
    absfile = '{0:s}.{1:s}.{2:s}.i{3:s}.abs_cells.dat'.format(galID,expn,ion,inc)
    gasfile = '{0:s}_GZa{1:s}.{2:s}.txt'.format(galID,expn,ion)
    
    # Read in the absorbing cells
    losnum, abscells = np.loadtxt(testLoc+absfile, skiprows=1, usecols=(0,2), 
                                unpack=True)


    # Read in the gas file
    xLoc, yLoc, zLoc, density, temperature, alphaZ = np.loadtxt(testLoc+gasfile,
                                                     skiprows=2, 
                                                     usecols=(1,2,3,7,8,15),
                                                     unpack=True)
    # Read in the details about the line of sight
    xen, yen, zen, xex, yex, zex = np.loadtxt(testLoc+'lines.dat', skiprows=2, 
                                    usecols=(0,1,2,3,4,5), unpack=True)
    # Read in impact parameters
    b = np.loadtxt(testLoc+'lines.info', skiprows=2, usecols=(1,), unpack=True)
    

    # Find the unique LOS numbers in abs cell file
    uniqLOS, uniqCounts = np.unique(losnum, return_counts=True)

    deviations = np.zeros(len(uniqLOS))

    # Get the distance along the LOS for each cell
    xs = []
    maxdist = []
    mindist = []
    spread = []
    for i, los in enumerate(uniqLOS):
        ind = int(los)-1

        s = []
        for j in range(len(abscells)):
            if losnum[j]==los: 
                cellid = int(abscells[j])
                # Get the cell properties
                x = xLoc[cellid-1]
                y = yLoc[cellid-1]
                z = zLoc[cellid-1]
        
                # Get the distance from the LOS entry point
                # to this cell
                s.append(get_delta_s(xen[ind], yen[ind], zen[ind], 
                                xex[ind], yex[ind], zex[ind], 
                                x, y, z, b[ind]) )
    
        ypoints = [los for point in s]
        ax3[ionnum].plot(s,ypoints, '.', color='k', ms=1)

        maxdist.append(max(s))
        mindist.append(min(s))
        spread.append(max(s)-min(s))
        mid = np.mean(s)
        for point in s:
            xs.append(point-mid)
            
#        dev = np.mean(s)
        # Get the standard deviation of this distribution
        dev = np.std(s)
        deviations[i] = dev

    density = kde.gaussian_kde(xs)
    xgrid = np.linspace( min(xs), max(xs) )
    ax2.plot(xgrid, density(xgrid), label=ion) 

    ax.hist(deviations, bins=numbins, range=(0,0.2), log = True, 
            histtype='step', label=ion)

    ax3[ionnum].set_xlabel('Distance along LOS')
    ax3[ionnum].set_ylabel('LOS Num')
    ax3[ionnum].set_ylim([0,1000])
    ax3[ionnum].set_xlim([0,1])
    ax3[ionnum].set_title(ion)

    print 'Max maxdist: {0:f}'.format(max(maxdist))
    print 'Min mindist: {0:f}'.format(min(mindist))
    print 'Min spread:  {0:f}'.format(min(spread))
    print 'Max spread:  {0:f}'.format(max(spread))
    
ax.set_xlabel('Standard Deviation of Cell Location Along LOS')
ax.set_ylabel('Counts')
ax.legend(frameon=False)
fig.savefig('deviations.pdf',bbox_inches='tight')

ax2.set_xlim([-0.1,0.1])
ax2.legend(frameon=False)
fig2.savefig('kde.pdf', bbox_inches='tight')

fig3.tight_layout()
fig3.savefig('absorbingCells.png', bbox_inches='tight')













