
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

for ion in ions:

    fout = open('{0:s}_meanStdInClusters.out'.format(ion), 'w')
    fout.write('k\tTemp\tDense\tSNII\n')
    print '\nIon = ',ion
    filename = '{0:s}_abscells.subset'.format(ion)
    cells = np.loadtxt(filename, skiprows=1)
    filename = filename.replace('abscells', 'abscells_full')
    x, y, z, dense, temp, snII = np.loadtxt(filename, skiprows=1, 
                                usecols=(0,1,2,3,4,5), unpack=True)

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
    maxClusterNum = 100
    for k in range(1,maxClusterNum+1):

        print 'Number of clusters = ',k
        # Retrieve the clusters
        cluster = fcluster(z, k, criterion='maxclust')
        meanks.append(k)
        for clusterNum in range(1,k+1):
            print '\tFor cluster number {0:d}, {1:d} members'.format(clusterNum,
                            np.sum(cluster[cluster==clusterNum])/2)
            groupN = dense[cluster==clusterNum]
            groupT = temp[cluster==clusterNum]
            groupZ = snII[cluster==clusterNum]
            
            ks.append(k)
            tempMean.append(np.log10(np.mean(groupT)))
            denseMean.append(np.log10(np.mean(groupN)))
            metalMean.append(np.log10(np.mean(groupZ)))

            tempSpread.append(np.log10(np.std(groupT)))
            denseSpread.append(np.log10(np.std(groupN)))
            metalSpread.append(np.log10(np.std(groupZ)))
            
            tempMad.append(np.log10(mad(groupT)))
            denseMad.append(np.log10(mad(groupN)))
            metalMad.append(np.log10(mad(groupZ)))

        fout.write('{0:d}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n'.format(k,
                    np.mean(tempSpread[-k]),
                    np.mean(denseSpread[-k]),
                    np.mean(metalSpread[-k]) ))
        meanStdT.append(np.mean(tempSpread[-k]))
        meanStdN.append(np.mean(denseSpread[-k]))
        meanStdZ.append(np.mean(metalSpread[-k]))

        meanMadT.append(np.mean(tempMad[-k]))
        meanMadN.append(np.mean(denseMad[-k]))
        meanMadZ.append(np.mean(metalMad[-k]))

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














