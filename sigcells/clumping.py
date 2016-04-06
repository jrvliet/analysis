
'''
A code to quantify the clumping of absorbing cells in spatial location, 
density, temperature, and metallicity
Makes files called std<property>.dat
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

expn = '0.510'
inc = '10'

numbins = 20

# Intiialize arrays    
stddevs = []
stddevs.append(ions)
for ion in ions:
    stddevs.append([])

spreads = []
spreads.append(ions)
for ion in ions:
    spreads.append([])

stdTemps = []
stdTemps.append(ions)
for ion in ions:
    stdTemps.append([])

stdDense = []
stdDense.append(ions)
for ion in ions:
    stdDense.append([])

stdMetal = []
stdMetal.append(ions)
for ion in ions:
    stdMetal.append([])

for galID, expn in zip(galIDs, expns):
    print ''
    print galID    

    for ionnum, ion in enumerate(ions):
        
        listfile = '{0:s}/{1:s}.{2:s}.list'.format(absLoc,galID,ion)
        listf = open(listfile, 'r')

        outfile = '{0:s}_{1:s}_abscellClumping.dat'.format(galID, ion)
        f = open(outfile, 'w')

        header = ('{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\t{5:s}\t'
                  '{6:s}\t{7:s}\t{8:s}\t{9:s}\t{10:s}\n')
        header = header.format('Expn','Min LOS length', 'Max LOS length', 
            'Mean LOS length', 'Min Dist', 'Max Dist', 'Min Spread', 
            'Max Spread', 'Mean Spread', 'Median Spread', 'Std Dev Spread')

        f.write(header)
        sFormat = ('{0:s}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{5:.6f}\t'
                   '{6:.6f}\t{7:.6f}\t{8:.6f}\t{9:.6f}\t{10:.6f}\n')

        print '\t\t',ion

        for line in listf:

            # Define files
            absfile = absLoc + line.strip()
            abssplit = absfile.split('.')
            a = '{0:s}.{1:s}'.format(abssplit[1],abssplit[2])
            gasfile = '{0:s}/{1:s}_GZa{2:s}.{3:s}.txt'.format(gasLoc,galID,a,ion)
            linesfile = '{0:s}/{1:s}_{2:s}_lines.dat'.format(linesLoc,galID,a)
            metalabsfile = '{0:s}/{1:s}_{2:s}_metalSelectedLOS.dat'.format(linesLoc,galID,a)

            # Read in the absorbing cells
            try:
                losnum, abscells = np.loadtxt(absfile, skiprows=1, 
                                          usecols=(0,2), unpack=True)
            except ValueError:
                continue

            # Read in the gas file
            try:
                xLoc, yLoc, zLoc, density, temperature, alphaZ = np.loadtxt(gasfile,
                                                             skiprows=2, 
                                                             usecols=(1,2,3,7,8,15),    
                                                         unpack=True)
            except ValueError:
                continue

            # Read in the details about the line of sight
            xen, yen, zen, xex, yex, zex = np.loadtxt(linesfile, 
                                skiprows=2, usecols=(0,1,2,3,4,5), unpack=True)
            # Read in impact parameters
#            b = np.loadtxt(loc+'lines.info', skiprows=2, usecols=(1,), unpack=True)
            

            # Find the unique LOS numbers in abs cell file
            uniqLOS, uniqCounts = np.unique(losnum, return_counts=True)

            deviations = np.zeros(len(uniqLOS))

            # Read in the LOS that have metal absoroption features detected
            metalLOS = np.loadtxt(metalabsfile, skiprows=1, usecols=(0,), unpack=True)

            #print 'Number of metal LOS: ',len(metalLOS)
#            print 'Number of unique Los: ',len(uniqLOS)
#            print 'Number of LOS in both arrys: ',len(np.intersect1d(metalLOS,uniqLOS))
#            print np.intersect1d(metalLOS, uniqLOS)
            #for i in range(0,len(metalLOS)):
            #    print metalLOS[i], uniqLOS[i]
#            sys.exit()
            # Get the distance along the LOS for each cell
            for i, los in enumerate(uniqLOS):
                if los in metalLOS:
                    ind = int(los)-1

                    s, l = [], []
                    ts, ns, zs = [], [], [],
                    for j in range(len(abscells)):
                        if losnum[j]==los: 
                            cellid = int(abscells[j])
                            # Get the cell properties
                            x = xLoc[cellid-1]
                            y = yLoc[cellid-1]
                            z = zLoc[cellid-1]
                            ts.append(temperature[cellid-1])
                            ns.append(density[cellid-1])
                            zs.append(alphaZ[cellid-1])
                    
                            # Get the distance from the LOS entry point
                            # to this cell
                            dist, leng = get_delta_s(xen[ind], yen[ind], zen[ind], 
                                            xex[ind], yex[ind], zex[ind], 
                                            x, y, z) 
                            s.append(dist)
                            l.append(leng)
                

                    # Get the standard deviation of this distribution
                    dev = np.std(s)
                    stddevs[ionnum+1].append(dev)
                    spreads[ionnum+1].append(max(s)-min(s))

                    stdTemps[ionnum+1].append(np.std(ts))
                    stdDense[ionnum+1].append(np.std(ns))
                    stdMetal[ionnum+1].append(np.std(zs))




#            result = sFormat.format(a,min(l),max(l),np.mean(l),
#                    max(maxdist),min(mindist),max(spread),min(spread),
#                    np.mean(spread),np.median(spread),np.std(spread))
#            f.write(result) 
        f.close()
 
stdfile = 'stdLocation.dat'
with open(stdfile, 'w') as f:
    json.dump(stddevs, f)


spreadfile = 'spreads.dat'
with open(spreadfile, 'w') as f:
    json.dump(spreads, f)

tempfile = 'stdTemperatures.dat'
with open(tempfile, 'w') as f:
    json.dump(stdTemps, f)

tempfile = 'stdDensities.dat'
with open(tempfile, 'w') as f:
    json.dump(stdDense, f)

tempfile = 'stdMetals.dat'
with open(tempfile, 'w') as f:
    json.dump(stdMetal, f)






