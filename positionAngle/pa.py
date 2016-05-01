
'''
Calculates the EW distribution as a funciton of 
position angle of the LOS around the galaxy.
'''

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats

def calc_pa(phi):
    '''
    Calculates the position angle for the LOS
    '''
    rot = math.floor(phi / 90)
    pa = phi - 90*rot
    return pa

def gmean(vals):
    '''
    Returns geometric mean of values in vals
    '''
    return stats.gmean(vals)


location = '/Users/jacob/research/velas/vela2b/vela23/a0.490/i90/'
absfile = 'vela2b-23.0.490.{0:s}.i90.abs_cells.dat'
sysabsfile = 'vela2b-23.{0:s}.a0.490.i90.ALL.sysabs'
ions = ['HI', 'MgII', 'CIV', 'OVI']
stat = 'mean' 
statname = 'mean'

b, angle = np.loadtxt(location+'lines.info', skiprows=2, usecols=(1,2), 
                        unpack=True)
phi = []
for p in angle:
    phi.append(calc_pa(p))

plt.hist(phi,bins=10,histtype='step')
plt.savefig('pa_hist.png', bbox_inches='tight')
plt.cla()
plt.clf()

fig1, ((ax11, ax12), (ax13, ax14)) = plt.subplots(2,2,figsize=(10,10))
axes1 = (ax11, ax12, ax13, ax14)
fig2, ((ax21, ax22), (ax23, ax24)) = plt.subplots(2,2,figsize=(10,10))
axes2 = (ax21, ax22, ax23, ax24)

minpa, maxpa = 0, 90
minew, maxew = -2, 2
nbins = 10

print('Ion\tMin\tMax\tMean\tMedian')
string = '{0:s}\t{1:.3f}\t{2:.3f}\t{3:.3f}\t{4:.3f}'
for ion, ax1, ax2 in zip(ions,axes1,axes2):

    ew,logN = np.loadtxt(location+sysabsfile.format(ion), skiprows=1, 
                    usecols=(5,19), unpack=True)
    pa, e, n = [], [], []
    for i in range(0,len(ew)):
        if ew[i]>0 and logN[i]>0:
            pa.append(phi[i])
            e.append(ew[i])
            n.append(logN[i])

    print(string.format(ion,np.min(logN),np.max(logN),np.mean(logN),
            np.median(logN)))

    bin_means, bin_edges, binnumber = stats.binned_statistic(pa, ew,
        statistic=stat, bins=nbins)
    ax1.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=2,
               label='binned statistic of data')
    ax1.vlines(bin_edges[1:-1], bin_means[:-1], bin_means[1:], colors='g')
    ax1.set_xlabel('PA [deg]')
    ax1.set_ylabel('EW [$\AA$]')
    ax1.set_title(ion)


    bin_means, bin_edges, binnumber = stats.binned_statistic(pa, logN,
        statistic=stat, bins=nbins)
    ax2.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=2,
               label='binned statistic of data')
    ax2.vlines(bin_edges[1:-1], bin_means[:-1], bin_means[1:], colors='g')
    ax2.set_xlabel('PA [deg]')
    ax2.set_ylabel('EW [$\AA$]')
    ax2.set_title(ion)


fig1.tight_layout()
fig1.savefig('pa_vela2b-23_ew_{0:s}.png'.format(statname), bbox_inches='tight')

fig2.tight_layout()
fig2.savefig('pa_vela2b-23_logN_{0:s}.png'.format(statname), bbox_inches='tight')





