
'''
A code to quantify the clumpin gof absorbing cells
based on their spatial location
'''

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from scipy.stats import kde
from jenks import jenks
import subprocess as sp
import json
import sys
import scipy.cluster.hierarchy as sh
import scipy.spatial.distance as sd
import time



def num_clusters(z):

    last = z[-10:,2]
    acceleration = np.diff(last, 2)  # 2nd derivative of the distances
    acceleration_rev = acceleration[::-1]
    k = acceleration_rev.argmax() + 2 

    return k

def mad(x):
    '''
    Computes the Mean Absolute Deviation of the sample
    Defined as:
        mad = median( |x - median(x)| )
    '''
    
    # Find the median of the sample
    med = np.median(x)
    
    # Get the absolute deviation
    dev = [abs(i-med) for i in x]
    
    # Compute the MAD
    val = np.median(dev)
    return val

def cov(x):
    '''
    Computes the coefficient of variation of the sample
    Defined as:
        cov = std(x)/mean(x)
    '''
    stdev = np.std(x)
    mean = np.mean(x)
    if stdev<0 or mean<0:
        print stdev, mean, min(x), max(x)
    coeff = np.std(x)/np.mean(x)
    return coeff

def group_radius(gx,gy,gz):
    '''
    Computes the radius of the group based on the
    x,y,z coordinates of the group members
    
    Represent the cluster as a sphere
    Diameter is mean of the largest distance between
    member points along each axis
    '''

    longestX = max(gx) - min(gx)
    longestY = max(gy) - min(gy)
    longestZ = max(gz) - min(gz)
    
    diameter = (longestX + longestY + longestZ) / 3.0

    return diameter / 2.0

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = sh.dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata





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


# Loop over galID
for galID, expn in zip(galIDs, expns):

    print galID
    
    # Loop over ions
    for ionnum, ion in enumerate(ions):

        # Open the files that contain the list of all abs_cells files
        listfile = '{0:s}/{1:s}.{2:s}.list'.format(absLoc,galID,ion)
        
        with open(listfile, 'r') as listf:
    
            for line in listf:

                # Open each abs_cell file and read it in
                absfile = line.strip()
                abssplit = absfile.split('.')
                a = '{0:s}.{1:s}'.format(abssplit[1],abssplit[2])
                gasfile = '{0:s}/{1:s}_GZa{2:s}.{3:s}.txt'.format(gasLoc,galID,a,ion)
                linesfile = '{0:s}/{1:s}_{2:s}_lines.dat'.format(linesLoc,galID,a)
    
                # Read in the absorbing cells
                try: 
                    abscells = np.loadtxt(absLoc+absfile, skiprows=1, usecols=(2,), unpack=True)
                except ValueError:
                    continue
    
                # Read in the gas file
                try:
                    xLoc, yLoc, zLoc, density, temperature, snII = np.loadtxt(gasfile,
                            skiprows=2, usecols=(1,2,3,7,8,9), unpack=True)
                except ValueError:
                    continue

                numCells = len(xLoc)
                numSubset = int(numCells/100.0)

                print 'Number of Cells = {0:d}\tnumSubset = {1:d}'.format(numCells, numSubset)

                # Select a number of random cells equal to the numSubset from the
                # the absorbing cells
                subsetLoc = np.zeros((numSubset, 3))
                subsetProp = np.zeros((numSubset, 3))
                
                for i in range(numSubset):
                    index = np.random.randint(0,numCells)
                    subsetLoc[i, 0] = xLoc[index]
                    subsetLoc[i, 1] = yLoc[index]
                    subsetLoc[i, 2] = zLoc[index]
    
                    subsetProp[i, 0] = density[index]
                    subsetProp[i, 1] = temperature[index]
                    subsetProp[i, 2] = snII[index]
                    

                # Genearate the linkage matrix
                # Use the Ward variance minimization algorithm
                time1 = time.time()
                z = sh.linkage(subsetLoc, 'ward')
                time2 = time.time()
                print 'Duration of linkage = {0:f} seconds'.format(time2-time1)

                # Determine how well the clustering preserves the
                # original distance
                c, coph_dist = sh.cophenet(z, sd.pdist(subsetLoc))
                print c


                # Create the fancy dendrogram and save
                fancy_dendrogram(z, truncate_mode='lastp', p=12,
                                leaf_rotation=90, leaf_font_size=12,
                                show_contracted=True, annotate_above=10)

                figname = '{0:s}_{1:s}_{2:s}_abscells_dendrogram.png'.format(ion,galID,expn)
                plt.savefig(figname, bbox_inches='tight')
                plt.cla()
                plt.clf()

                k = num_clusters(z)
                print 'Number of clusters = {0:d}'.format(k)



                sys.exit()
            
























