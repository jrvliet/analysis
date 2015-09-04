#!/usr/bin/python
# Filename: nT_raw_mass.py
#
# Author: Jacob Vander Vliet
# Version: 1
# Date: 4/12/2013
# Last Modified: 4/12/13
#
# Description:
#   Creates the phase diagram of the raw snapshots
#   The bins contain the total mass in each bin
#
# Instructions:
#   python nT_raw_mass.py
#

import numpy as np
import matplotlib.pyplot as plt
from os import system
import scipy.interpolate
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib as mpl
import sys

params = {'xtick.labelsize':11,
          'ytick.labelsize':11,
          'font.size':12}
plt.rcParams.update(params)

# Assume a constant mass fraction for hydrogen
Xh = 0.7

# Constants
Ah = 1.00794
AMU = 1.6605e-24  # grams 
MSUN = 1.99e33  # grams
PC = 3.086e18  # cm

galID_list = ['D9o2', 'D9q', 'D9m4a']
feedback_list = ['dwSN', 'dwALL_1', 'dwALL_8']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]


file_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/'
numbins = 50
normalize = 0
histos = []
xed = []
yed = []

# Create the subplots
fig,((p11,p12,p13)) = plt.subplots(1,3,figsize=(10.2, 3.5))
plot_list = [p11, p12, p13]

for i in range(0,len(galID_list)):

    galID = galID_list[i]
    expn = expn_list[i]

    n = []
    t = []
    m = []

    for a in expn:

        # Open the data. Cell size is in parsecs in the pre-rates boxes
        abs_file = galID+'_outputs/a'+a+'/'+galID+'_GZa'+a+'.txt'
        size, lognH, logT = np.loadtxt(file_loc+abs_file, skiprows=2, 
                                       usecols=(0, 7, 8), unpack=True)
        print galID, a

        for j in range(0,len(lognH)):

            n.append(np.log10(lognH[j]))
            t.append(np.log10(logT[j]))
            
            # Calculate the mass of each cell
            # Get the density from rho = Ah * AMU * density / Xh
            rho1 = Ah*AMU*lognH[j] / Xh      # grams per cubic cm
            
            # Convert to solar massses / cubic cm
            rho =  (rho1 / MSUN) / np.power(PC, -3.0)

            # Get the mass by multiplying by the cell volume
            volume = np.power(size[j], 3.0)
            mass = rho * volume
            m.append(mass)
            
            verbose = 0
            if verbose==1:
                print '\n Mass Calculation: '
                print '\t Density            : {0:e}'.format(lognH[j])
                print '\t AMU                : {0:e}'.format(AMU)
                print '\t Atomic Weight      : {0:e}'.format(Ah)
                print '\t Mass Fraction      : {0:e}'.format(Xh)
                print '\t Density [g/cm3]    : {0:e}'.format(rho1)
                print ''
                print '\t MSUN               : {0:e}'.format(MSUN)
                print '\t PC                 : {0:e}'.format(PC)
                print '\t Density [Msun/pc3] : {0:e}'.format(rho)
                print ''
                print '\t Volume [pc3]       : {0:e}'.format(volume)
                print '\t Mass               : {0:e}'.format(mass)

                sys.exit()

            
    H, xedges, yedges = np.histogram2d( n, t, bins=numbins, normed=normalize)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]


    dum_ind = 43
    print len(xedges)
    print len(yedges)
    print 'Index: ', xedges[dum_ind]
    print 'Histogram: '
    print H[dum_ind,:]

    # Change the value of each bin to contain the total mass in the cells
    # with that density and temperature instead of the number of cells
    print ''
    mass_arr = np.zeros_like(H)

    # Zero out the array
    for x_ind in range(0,numbins):
        for y_ind in range(0,numbins):
            mass_arr[x_ind, y_ind] = 0.




    xs = np.digitize(n, xedges)
    ys = np.digitize(t, yedges)


#    print [item for item in xs]
#    print max(xs)
#    print min(xs)


    for j in range(0,len(xs)):
        if xs[j]>49:
            xs[j]=49

    for j in range(0,len(ys)):
        if ys[j]>49:
            ys[j]=49


    for x, y, ms in zip(xs, ys, m):
        mass_arr[x-1, y-1] += ms


#    for j in range(0,len(m)):

        # Determine the index where this cell's density resides
#        density = n[j]
#        temperature = t[j]

#        n_index = np.digitize(density, xedges)
#        t_index = np.digitize(density, yedges)

#        for k in range(0,numbins):
#            if density>xedges[k] and density<xedges[k+1]:
#                n_index = k
#            if temperature>yedges[k] and temperature<yedges[k+1]:
#                t_index = k
        # Determine the index where this cell's temperature resides

#        for k in range(0,numbins):
#        dumflag = False
#        if H[n_index, t_index] == 0 and dumflag==True:
#            print 'ERROR\n\t density: {0:f} \t temp: {1:f} \t mass: {2:f} \t Index: ({3:d}, {4:d})'.format(density, temperature, m[j], n_index, t_index)
#        mass_arr[n_index, t_index] += m[j]
        
        

    print '\n Masses: '
    print mass_arr[dum_ind,:]

    # Rotate and masking
    mass_arr = np.rot90(mass_arr)
    mass_arr = np.flipud(mass_arr)
    Hmasked = np.ma.masked_where(mass_arr==0,mass_arr)
        
    # Take the log of the count
    Hmasked = np.log10(Hmasked)
    xed.append(xedges)
    yed.append(yedges)
    histos.append(Hmasked)

labels = ['(a)', '(b)', '(c)']
print ''
print len(histos)
print ''
for i in range(0,len(plot_list)):
    ax = plot_list[i]
    H = histos[i]
    xedges = xed[i]
    yedges = yed[i]
    mesh = ax.pcolormesh(xedges,yedges,H)
    ax.set_xlim([-8, 1])
    ax.set_ylim([2,8])
    ax.text(-1, 7, labels[i])

    ax.set_ylabel('$\log$ (T [K] )')
    ax.set_xlabel(' $\log (n_{H} $ [cm$^{-3}$] )')

    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.set_label('$\log$ (mass)')
#    plt.savefig('dum_'+str(i)+'.eps')
    
#dropx = [p12, p22, p32, p13, p23, p33, p11, p21, p31]
#plt.setp([a.get_xticklabels() for a in dropx],visible=False)

plt.tight_layout()
s = file_loc+'phase_bulk_masses.eps'
plt.savefig(s, bbox_inches='tight')
