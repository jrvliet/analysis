
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
    if phi<=90:
        return phi
    elif phi>90 and phi<=180:
        return 180.0 - phi
    elif phi>180 and phi<=270:
        return phi - 180.0
    else:
        return 350.0 - phi

def gmean(vals):
    '''
    Returns geometric mean of values in vals
    '''
    return stats.gmean(vals)


location = '/Users/jacob/research/velas/vela2b/vela23/a0.490/i90/'
location = '/home/jacob/research/velas/vela2b/vela28/a0.490/i90/'
absfile = 'vela2b-28.0.490.{0:s}.i90.abs_cells.dat'
sysabsfile = 'vela2b-28.{0:s}.a0.490.i90.ALL.sysabs'
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
fig3, ((ax31, ax32), (ax33, ax34)) = plt.subplots(2,2,figsize=(10,10))
axes3 = (ax31, ax32, ax33, ax34)
fig4, ((ax41, ax42), (ax43, ax44)) = plt.subplots(2,2,figsize=(10,10))
axes4 = (ax41, ax42, ax43, ax44)

minpa, maxpa = 0, 90
minew, maxew = -2, 2
nbins = 10

print('Ion\tMin\tMax\tMean\tMedian')
string = '{0:s}\t{1:.3f}\t{2:.3f}\t{3:.3f}\t{4:.3f}'
for ion, ax1, ax2, ax3, ax4 in zip(ions,axes1,axes2,axes3,axes4):

    sysfile = '{0:s}/{1:s}/{2:s}'.format(location,ion,sysabsfile.format(ion))
    ew,logN = np.loadtxt(sysfile, skiprows=1,
                    usecols=(5,19), unpack=True)
    pa, e, n = [], [], []
    for i in range(0,len(ew)):
        if ew[i]>0 and logN[i]>0:
            pa.append(phi[i])
            e.append(ew[i])
            n.append(logN[i])

    print(string.format(ion,np.min(ew),np.max(ew),np.mean(ew),
            np.median(ew)))

    bin_means, bin_edges, binnumber = stats.binned_statistic(pa, e,
        statistic=stat, bins=nbins)
    ax1.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=2,
               label='binned statistic of data')
    ax1.vlines(bin_edges[1:-1], bin_means[:-1], bin_means[1:], colors='g')
    ax1.set_xlabel('PA [deg]')
    ax1.set_ylabel('EW [$\AA$]')
    ax1.set_title(ion)


    bin_means, bin_edges, binnumber = stats.binned_statistic(pa, n,
        statistic=stat, bins=nbins)
    ax2.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=2,
               label='binned statistic of data')
    ax2.vlines(bin_edges[1:-1], bin_means[:-1], bin_means[1:], colors='g')
    ax2.set_xlabel('PA [deg]')
    ax2.set_ylabel('EW [$\AA$]')
    ax2.set_title(ion)


    # Make a polar plot 
    # being mean of EW at that location
    x, y, equiv = [], [], []
    for r, theta, width in zip(b,np.radians(phi),ew):
        x.append(r*np.cos(theta))
        y.append(r*np.sin(theta))
        equiv.append(width)
    ax3.scatter(x, y, c=ew, s=10, marker='o')
    ax3.set_xlim([0.0,r.max()])
    ax3.set_ylim([0.0,r.max()])
    ax3.set_title(ion)

    # Make 2d histogram of x, y position on sky, with bin value
    bin_means, x_edges, y_edges, binnumber = stats.binned_statistic_2d(x,y,equiv,
            statistic='mean', bins=50)
    whereNans = np.isnan(bin_means)
    bin_means[whereNans] = 0.0
    print ew.min(), ew.max()
    print bin_means.min(), bin_means.max()
    bin_means = np.ma.masked_where(bin_means==0, bin_means)
    bin_means = np.log10(bin_means)
    mesh = ax4.pcolormesh(x_edges, y_edges, bin_means)
    ax4.set_xlim([0,x_edges[-1]])
    ax4.set_ylim([0,y_edges[-1]])
    ax4.set_title(ion)
    cbar = fig4.colorbar(mesh, ax=ax4, use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('Mean EW', rotation=270, fontsize=12)

fig1.tight_layout()
fig1.savefig('pa_vela2b-23_ew_{0:s}.png'.format(statname), bbox_inches='tight')

fig2.tight_layout()
fig2.savefig('pa_vela2b-23_logN_{0:s}.png'.format(statname), bbox_inches='tight')

fig3.tight_layout()
fig3.savefig('polar_vela2b-23_ew_{0:s}.png'.format(statname), bbox_inches='tight')

fig4.tight_layout()
fig4.savefig('pahist_vela2b-23_ew_{0:s}.png'.format(statname), bbox_inches='tight')





