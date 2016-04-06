
'''
A code to determine the spread in LOS location vs LOS velocit
Writes <ion>_vlos.dat
'''
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from scipy.stats import kde
from jenks import jenks
import subprocess as sp
import json
import sys

def get_delta_s(xen, yen, zen, xex, yex, zex, x, y, z):
    
    '''
    Returns the distance the point given by (x,y,z) is from
    (xen, yen, zen) along the LOS defined by the starting
    and ending points
    '''
    
    # Length of LOS
    l = sqrt( (xen-xex)**2 + (yen-yex)**2 + (zen-zex)**2 )

    # Distance from entry point to cell
    s = sqrt( (xen-x)**2 + (yen-y)**2 + (zen-z)**2 )

    return s, l

def get_vlos(xen, yen, zen, xex, yex, zex, vx, vy, vz, l):
    '''
    Returns the component of the velocity vector along the 
    direction of the LOS, as found by the dot product of 
    the directional vector and the velocity vector
    '''

    # Normalized direcitonal vector
    dx = (xex-xen) / l
    dy = (yex-yen) / l
    dz = (zex-zen) / l

    # Find the dot product
    vlos = dx*vx + dy*vy + dz*vz
    
    return vlos



testLoc = '/home/hyades/jrvander/exampleData/'
baseLoc = '/home/matrix3/jrvander/sebass_gals/dwarfs/'
absLoc = '/home/jacob/research/dwarfs/abscells/individual/'
gasLoc = '/home/jacob/research/dwarfs/gasfiles/'
linesLoc = '/home/jacob/research/dwarfs/lines/'

ions = ['HI', 'MgII', 'CIV', 'OVI']

galIDs = ['D9o2', 'D9q', 'D9m4a']
expn1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']

expn1 = ['1.002']
expn2 = ['1.001']
expn3 = ['1.000']
expns = [expn1, expn2, expn3]

for ionnum, ion in enumerate(ions):
    print ion
    s, l, v = [], [], []
    f = open('{0:s}_vlos.dat'.format(ion), 'w')
    header = 'LOS length\tS\t\tVlos\t\tDensity\t\tTemp\t\tAlphaZ\t\tCellSize\tnIon\n'
    f.write(header)
    form = '{0:.6f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{5:.6f}\t{6:.6f}\t{7:.6f}\n'
    for galID in galIDs:
        print '\t',galID    

        listfile = '{0:s}/{1:s}.{2:s}.list'.format(absLoc,galID,ion)
        listf = open(listfile, 'r')

        for line in listf:
    
            # Set file names
            absfile = absLoc + line.strip()
            abssplit = absfile.split('.') 
            a = '{0:s}.{1:s}'.format(abssplit[1],abssplit[2])
            gasfile = '{0:s}/{1:s}_GZa{2:s}.{3:s}.txt'.format(gasLoc,galID,a,ion)
            linesfile = '{0:s}/{1:s}_{2:s}_lines.dat'.format(linesLoc,galID,a)

            # Read in the absorbing cells
            try:
                losnum, abscells = np.loadtxt(absfile, skiprows=1, 
                                          usecols=(0,2), unpack=True)
            except ValueError:
                continue

            # Read in the gas file
            try:
                size, xLoc, yLoc, zLoc, xVel, yVel, zVel, density, temperature, nIon, alphaZ = np.loadtxt(
                                                             gasfile,
                                                             skiprows=2, 
                                                             usecols=(0,1,2,3,4,5,6,7,8,13,15),    
                                                             unpack=True)
            except ValueError:
                continue

            # Read in the details about the line of sight
            xen, yen, zen, xex, yex, zex = np.loadtxt(linesfile, 
                                skiprows=2, usecols=(0,1,2,3,4,5), unpack=True)
            
            # Find the unique LOS numbers in abs cell file
            uniqLOS, uniqCounts = np.unique(losnum, return_counts=True)

            deviations = np.zeros(len(uniqLOS))

            # Get the distance along the LOS for each cell
            for i, los in enumerate(uniqLOS):
                ind = int(los)-1

                for j in range(len(abscells)):
                    if losnum[j]==los: 
                        cellid = int(abscells[j])
                        # Get the cell properties
                        x = xLoc[cellid-1]
                        y = yLoc[cellid-1]
                        z = zLoc[cellid-1]
                        vx = xVel[cellid-1]
                        vy = yVel[cellid-1]
                        vz = zVel[cellid-1]
                        ionDense = np.log10(nIon[cellid-1])
                        n = np.log10(density[cellid-1])
                        t = np.log10(temperature[cellid-1])
                        metal = alphaZ[cellid-1]
                        cellL = size[cellid-1]
                
                        # Get the distance from the LOS entry point
                        # to this cell
                        dist, leng = get_delta_s(xen[ind], yen[ind], zen[ind], 
                                        xex[ind], yex[ind], zex[ind], 
                                        x, y, z) 
                        vlos = get_vlos(xen[ind], yen[ind], zen[ind],
                                        xex[ind], yex[ind], zex[ind],
                                        vx, vy, vz, leng)

                        s = dist/leng
                        f.write(form.format(leng, s, vlos, n, t, metal, cellL, ionDense))
        listf.close()
    f.close()





















