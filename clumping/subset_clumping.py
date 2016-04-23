
'''
Code to determine the clumping of absorbing cells 
Uses only the subset files: <ion>_abscells.subset
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster
from scipy.spatial.distance import pdist
import sys

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

    ddata = dendrogram(*args, **kwargs)

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

ions = ['HI', 'MgII', 'CIV', 'OVI']
plotting = 0

for ion in ions:

    fout = open('{0:s}_meanStdInClusters.out'.format(ion), 'w')
    fout.write('k\tTemp\tDense\tSNII\n')
    print '\nIon = ',ion
    filename = '{0:s}_abscells.subset'.format(ion)
    cells = np.loadtxt(filename, skiprows=1)
    filename = filename.replace('abscells', 'abscells_full')
    xloc, yloc, zloc, dense, temp, snII = np.loadtxt(filename, skiprows=1, 
                                usecols=(0,1,2,3,4,5), unpack=True)

    print '\nLimits: '
    print '\tDensity:     {0:.4e}\t{1:.4e}'.format(min(dense),max(dense))
    print '\tTemperature: {0:.4e}\t{1:.4e}'.format(min(temp),max(temp))
    print '\tMetallicity: {0:.4e}\t{1:.4e}'.format(min(snII),max(snII))
    # Generate the linkage matrix
    # Use the Ward variance minimization algorithm
    z = linkage(cells, 'ward')

    # Determine the how well the clustering preserves the 
    # original distance
    c, coph_dist = cophenet(z, pdist(cells))
    print c

    # Save the linking matrix to file
    zfilename = '{0:s}_abscells.linking'.format(ion)
    header = 'Index 1\tIndex 2\tDistance\tNum in Cluster\n'
    np.savetxt(zfilename, z, delimiter='\t', header=header)


    # Create the dendrogram and save
    fancy_dendrogram(z, truncate_mode='lastp', p=12,
                    leaf_rotation=90, leaf_font_size=12,
                    show_contracted=True, annotate_above=10)
    figname = '{0:s}_abscells_dendrogram.png'.format(ion)
    plt.savefig(figname, bbox_inches='tight')
    plt.cla()
    plt.clf()


    last = z[-10:, 2]
    last_rev = last[::-1]
    idxs = np.arange(1, len(last) + 1)
    plt.plot(idxs, last_rev)

    acceleration = np.diff(last, 2)  # 2nd derivative of the distances
    acceleration_rev = acceleration[::-1]
    plt.plot(idxs[:-2] + 1, acceleration_rev)
    figname = '{0:s}_abscells_accel.png'.format(ion)
    plt.savefig(figname, bbox_inches='tight')
    plt.cla()
    plt.clf()
    k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters
    print "clusters:", k

    # Loop over the number of clusters
    ks, meanks = [], []
    tempSpread, denseSpread, metalSpread = [], [], []
    tempMean, denseMean, metalMean = [], [], []
    tempMad, denseMad, metalMad = [], [], [] 
    meanStdT, meanStdN, meanStdZ = [], [], []
    meanMadT, meanMadN, meanMadZ = [], [], []
    maxClusterNum = 10

    ftempname = '{0:s}_clump_temperature_cov.out'.format(ion)
    ftemp = open(ftempname, 'w')    

    fnumberName = '{0:s}_clump_member_count.out'.format(ion)
    fnumber = open(fnumberName, 'w')

    fdenseName = '{0:s}_clump_density_cov.out'.format(ion)
    fdense = open(fdenseName, 'w')

    fmetalName = '{0:s}_clump_snII_cov.out'.format(ion)
    fmetal = open(fmetalName, 'w')

    fradiusName = '{0:s}_clump_radius.out'.format(ion)
    fradius = open(fradiusName, 'w')

    for k in range(1,maxClusterNum+1):

        tempstr = '{0:d}\t'.format(k)
        numstr = '{0:d}\t'.format(k)
        densestr = '{0:d}\t'.format(k)
        metalstr = '{0:d}\t'.format(k)
        radiusstr = '{0:d}\t'.format(k)

        # Retrieve the clusters
        cluster = fcluster(z, k, criterion='maxclust')
        meanks.append(k)
        for clusterNum in range(1,k+1):

            numberOfMembers = len(cluster[cluster==clusterNum])
            #print '\tFor cluster number {0:d}, {1:d} members'.format(clusterNum,
            #            numberOfMembers)
            groupN = dense[cluster==clusterNum]
            groupT = temp[cluster==clusterNum]
            groupM = snII[cluster==clusterNum]
            groupX = xloc[cluster==clusterNum]
            groupY = yloc[cluster==clusterNum]
            groupZ = zloc[cluster==clusterNum]
            
            groupRadius = group_radius(groupX, groupY, groupZ)
            ks.append(k)

            numstr += '{0:d}\t'.format(numberOfMembers)
            tempstr  += '{0:.4e}\t'.format(cov(groupT))
            densestr += '{0:.4e}\t'.format(cov(groupN))
            metalstr += '{0:.4e}\t'.format(cov(groupM))
            radiusstr += '{0:.4e}\t'.format(groupRadius)

        fnumber.write(numstr+'\n')
        ftemp.write(tempstr+'\n')
        fdense.write(densestr+'\n')
        fmetal.write(metalstr+'\n')
        fradius.write(radiusstr+'\n')

#        fout.write('{0:d}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n'.format(k,
#                    np.mean(tempSpread[-k]),
#                    np.mean(denseSpread[-k]),
#                    np.mean(metalSpread[-k]) ))
    ftemp.close()
    fnumber.close()
    fdense.close()
    fmetal.close()
    fradius.close()

    if plotting==1:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
        ax1.plot(ks, tempSpread, 'x')
        ax1.set_xlabel('Number of Clusters')
        ax1.set_ylabel('Std Dev of Temp in Each Cluster')
        ax1.set_xlim([0,maxClusterNum+1])

        ax2.plot(ks, denseSpread, 'x')
        ax2.set_xlabel('Number of Clusters')
        ax2.set_ylabel('Std Dev of Density in Each Cluster')
        ax2.set_xlim([0,maxClusterNum+1])

        ax3.plot(ks, metalSpread, 'x')
        ax3.set_xlabel('Number of Clusters')
        ax3.set_ylabel('Std Dev of SNII MF in Each Cluster')
        ax3.set_xlim([0,maxClusterNum+1])
        
        fig.tight_layout()
        fig.savefig('{0:s}_cloud_props_std.png'.format(ion), bbox_inches='tight')

        plt.cla()
        plt.clf()   
        plt.close('all')

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
        ax1.plot(ks, tempMean, 'x')
        ax1.set_xlabel('Number of Clusters')
        ax1.set_ylabel('Mean of Temp in Each Cluster')
        ax1.set_xlim([0,maxClusterNum+1])

        ax2.plot(ks, denseMean, 'x')
        ax2.set_xlabel('Number of Clusters')
        ax2.set_ylabel('Mean of Density in Each Cluster')
        ax2.set_xlim([0,maxClusterNum+1])

        ax3.plot(ks, metalMean, 'x')
        ax3.set_xlabel('Number of Clusters')
        ax3.set_ylabel('Mean of SNII MF in Each Cluster')
        ax3.set_xlim([0,maxClusterNum+1])
        
        fig.tight_layout()
        fig.savefig('{0:s}_cloud_props_mean.png'.format(ion), bbox_inches='tight')
        plt.cla()
        plt.clf()   
        plt.close('all')
            

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
        ax1.plot(ks, tempMad, 'x')
        ax1.set_xlabel('Number of Clusters')
        ax1.set_ylabel('MAD of Temp in Each Cluster')
        ax1.set_xlim([0,maxClusterNum+1])

        ax2.plot(ks, denseMad, 'x')
        ax2.set_xlabel('Number of Clusters')
        ax2.set_ylabel('MAD of Density in Each Cluster')
        ax2.set_xlim([0,maxClusterNum+1])

        ax3.plot(ks, metalMad, 'x')
        ax3.set_xlabel('Number of Clusters')
        ax3.set_ylabel('MAD of SNII MF in Each Cluster')
        ax3.set_xlim([0,maxClusterNum+1])
        
        fig.tight_layout()
        fig.savefig('{0:s}_cloud_props_mad.png'.format(ion), bbox_inches='tight')
        plt.cla()
        plt.clf()   
        plt.close('all')


        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
        ax1.plot(meanks, meanMadT, 'x')
        ax1.set_xlabel('Number of Clusters')
        ax1.set_ylabel('Mean MAD of Temp in Each Cluster')
        ax1.set_xlim([0,maxClusterNum+1])

        ax2.plot(meanks, meanMadN, 'x')
        ax2.set_xlabel('Number of Clusters')
        ax2.set_ylabel('Mean MAD of Density in Each Cluster')
        ax2.set_xlim([0,maxClusterNum+1])

        ax3.plot(meanks, meanMadZ, 'x')
        ax3.set_xlabel('Number of Clusters')
        ax3.set_ylabel('Mean MAD of SNII MF in Each Cluster')
        ax3.set_xlim([0,maxClusterNum+1])
        
        fig.tight_layout()
        fig.savefig('{0:s}_numClusterEffect_mad.png'.format(ion), bbox_inches='tight')
        plt.cla()
        plt.clf()   
        plt.close('all')
        

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
        ax1.plot(meanks, meanStdT, 'x')
        ax1.set_xlabel('Number of Clusters')
        ax1.set_ylabel('Mean Std of Temp in Each Cluster')
        ax1.set_xlim([0,maxClusterNum+1])

        ax2.plot(meanks, meanStdN, 'x')
        ax2.set_xlabel('Number of Clusters')
        ax2.set_ylabel('Mean Std of Density in Each Cluster')
        ax2.set_xlim([0,maxClusterNum+1])

        ax3.plot(meanks, meanStdZ, 'x')
        ax3.set_xlabel('Number of Clusters')
        ax3.set_ylabel('Mean Std of SNII MF in Each Cluster')
        ax3.set_xlim([0,maxClusterNum+1])
        
        fig.tight_layout()
        fig.savefig('{0:s}_numClusterEffect_std.png'.format(ion), bbox_inches='tight')
        plt.cla()
        plt.clf()   
        plt.close('all')
    
    fout.close()














