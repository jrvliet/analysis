
# Plots the EW distribution 

# Usage:
#   python 


import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import os
import sys
from scipy.optimize import curve_fit

def schechter(l, phi, lstar, alpha):

    schec = (phi/lstar) * pow((l/lstar), alpha) * np.exp(-1.0*l/lstar)
    return schec    

def fit_schech(xdata,ydata,errdata,index):
    
    p0 = [0.1, 2.0, -1.0]
    x = xdata[index:]
    y = ydata[index:]
    err = errdata[index:]
    
    p, covar = curve_fit(schechter, x, y, p0=p0, 
                         sigma=err, absolute_sigma=True,
                         bounds=(-np.inf, [np.inf, np.inf, 0.]),
                         method='dogbox')

    return p, covar

ion_list = ['HI', 'MgII', 'CIV', 'OVI']

# Read in the galaxy properties from galaxy.props
with open('galaxy.props','r') as f:
    galID = f.readline().split()[1]
    expn = f.readline().split()[1]
    redshift = f.readline().split()[1]
    mvir = f.readline().split()[1]
    rvir = float(f.readline().split()[1])
    inc = int(float(f.readline().split()[1]))

# Open output file
outfile = 'ewdist_schechter_{0:s}_a{1:s}_i{2:d}.out'.format(galID, expn, inc)
fout = open(outfile, 'w')
header = 'Ion\tPhi\t\t\tW*\t\t\tAlpha\n'
fout.write(header)

# Command line arguements for binfreq
column = 6              # Column the data are in (count from 1)
linear = 0              # Input data type (0=linear, 1=log10)
binsize = 0.1           # Binsize (equal log10)
lowerlimit = -3         # Lower bin limit (log10)
upperlimit = 2          # Upper bin limit (log10)
header = 1              # Number of header rows

# Get the location of the code
pathname = os.path.dirname(sys.argv[0])
codeLoc = os.path.abspath(pathname)

# Loop over ions
for ion in ion_list:
    
    allfile = './{0:s}/{1:s}.{0:s}.a{2:s}.i{3:d}.ALL.sysabs'.format(ion,galID,expn,inc)
    
    print ion
    # Run Chris's binning program
    command = codeLoc + '/bindata-logfreq {0:s} {1:d} {2:d} {3:f} {4:f} {5:f} {6:d}'.format(allfile, column, linear, binsize, lowerlimit, upperlimit, header)
    sp.call(command, shell=True)

    # Output will be named <galID>.logfreqbin
    # Rename to <galID>.<expn>.<ion>.ew.logfreqbin
    oldname = '.logfreqbin'.format(galID)
    newname = '{0:s}.{1:s}.{2:s}.ew.logfreqbin'.format(galID, expn, ion)
    command = 'mv {0:s} {1:s}'.format(oldname, newname)
    sp.call(command, shell=True)

    # Read in the binned data
    data = np.loadtxt(newname, skiprows=4)
    binCenter = data[:,0]
    freq = data[:,1]
    halfbin = data[:,2]
    errDown = -1.0*data[:,3]
    errUp   = data[:,4]

    # Fit a Schecter function to the data
    xdata, ydata, err = [], [], []
    unlogx, unlogy, unlogerr = [], [], []
    for i in range(0,len(freq)):
        if freq[i]>-50:
            xdata.append(binCenter[i]) 
            ydata.append(freq[i]) 
            err.append( (errUp[i]+errDown[i])/2.0 )
            unlogx.append(np.power(10.0, binCenter[i]))
            unlogy.append(np.power(10.0, freq[i]))
            unlogerr.append(np.power(10.0, (errUp[i]+errDown[i])/2.0))

    # Find the max of the distribution. There is a possibility of a turnover
    # at weak EW that should be ignored. This turnover is due to sensitivity
    # limits and is not indicitive of the actual results
    maxIndex = ydata.index(max(ydata))

    p, covar = fit_schech(unlogx, unlogy, unlogerr, maxIndex)
    print p
    sigma = np.sqrt(np.diag(covar))
    phi = p[0]
    lstar = p[1]
    alpha = p[2]

    # Write results to file
    s = '{0:s}\t{1:f} +/- {2:f}\t{3:f} +/- {4:f}\t{5:f} +/- {6:f}\n'.format(
        ion, phi, sigma[0], lstar, sigma[1], alpha, sigma[2])
    fout.write(s)

    # Plot the data
    subplotnum = 221+ion_list.index(ion)
    plt.subplot(subplotnum)
    xerrbin = pow(10.0,halfbin)
    yerrbinDown = pow(10.0, errDown)
    yerrbinUp = pow(10.0, errUp)
    plt.errorbar(binCenter, freq, xerr=halfbin, yerr=[errDown,errUp], 
                linestyle='none', label='Data')

    # Overplot fit
    x = []
    y = []
    for l in xdata:
        val = schechter(10**l, phi, lstar, alpha)
        y.append(np.log10(val))
        x.append(l)

    plt.plot(x, y, 'r', label='Fit')
    plt.xlabel('log( EW [$\AA$] )')
    plt.ylabel('log ( n(EW) )')
    plt.xlim([-3, 2])
    plt.ylim([-5, 2])

    plt.legend(frameon=False, loc='lower left', prop={'size':8})

plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.suptitle(r'{0:s}, a={1:s}, Rvir={2:.1f} kpc, i={3:d}'.format(galID, expn, rvir, inc))
s = '{0:s}_a{1:s}_i{2:d}_ewdist.pdf'.format(galID, expn, inc)
plt.savefig(s, bbox_inches='tight')
fout.close()





