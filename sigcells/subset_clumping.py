
'''
Code to determine the clumping of absorbing cells 
Uses only the subset files: <ion>_abscells.subset
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist


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

    print '\nIon = ',ion
    filename = '{0:s}_abscells.subset'.format(ion)
    cells = np.loadtxt(filename, skiprows=1)

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
