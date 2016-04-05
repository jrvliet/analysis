
'''
A code to quantify the clumping of absorbing cells in spatial location, 
density, temperature, and metallicity
'''
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from scipy.stats import kde
from jenks import jenks
import subprocess as sp
import json
import sys

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

    return s, l

testLoc = '/home/hyades/jrvander/exampleData/'
baseLoc = '/home/matrix3/jrvander/sebass_gals/dwarfs/'

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

for ion in ions:
    filename = '{0:s}_stddev.dat'.format(ion)
    command = 'touch {0:s}'.format(filename)
    sp.call(command, shell=True)
    command = command.replace('touch', 'rm')
    sp.call(command, shell=True)
    
stddevs = []
stddevs.append(ions)
for ion in ions:
    stddevs.append([])

spreads = []
spreads.append(ions)
for ion in ions:
    spreads.append([])

for galID, expn in zip(galIDs, expns):
    print ''
    print galID    

    for ionnum, ion in enumerate(ions):
        
        outfile = '{0:s}_{1:s}_abscellClumping.dat'.format(galID, ion)
        f = open(outfile, 'w')

#        header = ('{0:15>s}\t{1:15>s}\t{2:15>s}\t{3:15>s}\t{4:15>s}\t{5:15>s}\t'
#                  '{6:15>s}\t{7:15>s}\t{8:15>s}\t{9:15>s}\t{10:15>s}\n')
        header = ('{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\t{5:s}\t'
                  '{6:s}\t{7:s}\t{8:s}\t{9:s}\t{10:s}\n')
        header = header.format('Expn','Min LOS length', 'Max LOS length', 
            'Mean LOS length', 'Min Dist', 'Max Dist', 'Min Spread', 
            'Max Spread', 'Mean Spread', 'Median Spread', 'Std Dev Spread')

        f.write(header)
        sFormat = ('{0:s}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{5:.6f}\t'
                   '{6:.6f}\t{7:.6f}\t{8:.6f}\t{9:.6f}\t{10:.6f}\n')

        print '\t\t',ion
        for a in expn:
            print '\t', a

#            fig, ax = plt.subplots()

#            fig2, ((ax21, ax22),(ax23,ax24)) = plt.subplots(2,2)
#            ax2 = (ax21, ax22, ax23, ax24)

#            fig3, ((ax31, ax32),(ax33,ax34)) = plt.subplots(2,2)
#            ax3 = (ax31, ax32, ax33, ax34)
            loc = '{0:s}/{1:s}_outputs/a{2:s}/{3:s}/'.format(baseLoc, 
                                                             galID, a, ion)

            absfile = '{0:s}.{1:s}.{2:s}.abs_cells.dat'.format(galID,a,ion)
            gasfile = '{0:s}_GZa{1:s}.{2:s}.txt'.format(galID,a,ion)
    
            # Read in the absorbing cells
            try:
                losnum, abscells = np.loadtxt(loc+absfile, skiprows=1, 
                                          usecols=(0,2), unpack=True)
            except ValueError:
                continue
            # Read in the gas file
            try:
                xLoc, yLoc, zLoc, density, temperature, alphaZ = np.loadtxt(loc+gasfile,
                                                             skiprows=2, 
                                                             usecols=(1,2,3,7,8,15),    
                                                         unpack=True)
            except ValueError:
                continue
            # Read in the details about the line of sight
            xen, yen, zen, xex, yex, zex = np.loadtxt(loc+'../lines.dat', 
                                skiprows=2, usecols=(0,1,2,3,4,5), unpack=True)
            # Read in impact parameters
            b = np.loadtxt(loc+'lines.info', skiprows=2, usecols=(1,), unpack=True)
            

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

                s, l = [], []
                for j in range(len(abscells)):
                    if losnum[j]==los: 
                        cellid = int(abscells[j])
                        # Get the cell properties
                        x = xLoc[cellid-1]
                        y = yLoc[cellid-1]
                        z = zLoc[cellid-1]
                
                        # Get the distance from the LOS entry point
                        # to this cell
                        dist, leng = get_delta_s(xen[ind], yen[ind], zen[ind], 
                                        xex[ind], yex[ind], zex[ind], 
                                        x, y, z, b[ind]) 
                        s.append(dist/leng)
                        l.append(leng)
            
                ypoints = [los for point in s]
#                ax3[ionnum].plot(s,ypoints, '.', color='k', ms=1)
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
                stddevs[ionnum+1].append(dev)
                spreads[ionnum+1].append(max(s)-min(s))
            density = kde.gaussian_kde(xs)

            # Plot histogram of deviations from mean
#            ax.hist(deviations, bins=numbins, range=(0,0.2), log = True, 
#                    histtype='step', label=ion)
#            fig.suptitle('{0:s}, a{1:s}'.format(galID, a)
    
            # Plot histogram of spread
#            ax2[ionnum].hist(spread, bins=numbins, histtype='step', label=ion)
#            ax2[ionnum].set_xlabel('Spread')
#            ax2[ionnum].set_ylabel('Count')
#            ax2[ionnum].set_xlim([0,1])
#            ax2[ionnum].set_title(ion)
#            fig2.suptitle('{0:s}, a{1:s}'.format(galID, a)


            # Plot location of absorbing cells along los
#            ax3[ionnum].set_xlabel('Distance along LOS')
#            ax3[ionnum].set_ylabel('LOS Num')
#            ax3[ionnum].set_ylim([0,1000])
#            ax3[ionnum].set_xlim([0,1])
#            ax3[ionnum].set_title(ion)
#            fig3.suptitle('{0:s}, a{1:s}'.format(galID, a)

            result = sFormat.format(a,min(l),max(l),np.mean(l),
                    max(maxdist),min(mindist),max(spread),min(spread),
                    np.mean(spread),np.median(spread),np.std(spread))
            f.write(result) 
#            print 'Max maxdist:    {0:f}'.format(max(maxdist))
#            print 'Min mindist:    {0:f}'.format(min(mindist))
#            print 'Min spread:     {0:f}'.format(min(spread))
#            print 'Max spread:     {0:f}'.format(max(spread))
#            print 'Mean spread:    {0:f}'.format(np.mean(spread))
#            print 'Median spread:  {0:f}'.format(np.median(spread))
#            print 'Std Dev spread: {0:f}'.format(np.std(spread))

        f.close()
 
#    ax.set_xlabel('Standard Deviation of Cell Location Along LOS')
#    ax.set_ylabel('Counts')
#    ax.legend(frameon=False)
#    fig.savefig('deviations.pdf',bbox_inches='tight')


#    fig2.tight_layout()
#    fig2.savefig('kde.pdf', bbox_inches='tight')

#    fig3.tight_layout()
#    fig3.savefig('absorbingCells.png', bbox_inches='tight')



stdfile = 'stddev.dat'

with open(stdfile, 'w') as f:
    json.dump(stddevs, f)


spreadfile = 'spreads.dat'
with open(spreadfile, 'w') as f:
    json.dump(spreads, f)






sys.exit()

stdf = open(stdfile, 'w')
header = ''
for ion in ions:
    header += '{0:s}\t'.format(ion)
stdf.write(header+'\n')

# Write devs to file
import pdb; pdb.set_trace()
for j in range(0,len(stddevs[0])):
    s = ''
    for i in range(0,len(ions)):
        s += '{0:.6f}\t'.format(stddevs[i][j])
    
    stdf.write(s+'\n')

stdf.close()






