
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


ion_list = ['HI', 'MgII', 'CIV', 'OVI']
# Initial guess for Schechter parameters
# Phi = 1
# Lstar = -0.100
# Alpha = -1.0
p0 = [1.000, 1.000, -1.0]

# Read in the galaxy properties from galaxy.props
f = open('galaxy.props')
galID = f.readline().split()[1]
expn = f.readline().split()[1]
redshift = f.readline().split()[1]
mvir = f.readline().split()[1]
rvir = float(f.readline().split()[1])
f.close()


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
    
    allfile = './'+ion+'/'+galID+'.'+ion+'.a'+expn+'.ALL.sysabs'
    
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
    xdata = pow(10.0,binCenter)
    ydata = pow(10.0,freq)

    print p0 
    for i in range(0,len(xdata)):
        val = schechter(xdata[i], p0[0], p0[1], p0[2])
    (phi, lstar, alpha), paramCovariance = curve_fit(schechter, xdata, ydata, p0=p0)

    print 'Phi:\t{0:f}\nL*:\t{1:f}\nAlpha:\t{2:f}\n'.format(phi, lstar, alpha)
    print paramCovariance
    print 'One sigma :', np.sqrt(np.diag(paramCovariance))
    print ''

    # Plot the data
    subplotnum = 221+ion_list.index(ion)
    plt.subplot(subplotnum)
    xerrbin = pow(10.0,halfbin)
    yerrbinDown = pow(10.0, errDown)
    yerrbinUp = pow(10.0, errUp)
    plt.errorbar(xdata, ydata, xerr=halfbin, yerr=[errDown,errUp], 
                linestyle='none', label='Data')
#    plt.errorbar(binCenter, freq, xerr=halfbin, yerr=[errDown,errUp], 
#                linestyle='none', label='Data')


    # Overplot fit
    y = []
    for l in xdata:
        val = schechter(l, phi, lstar, alpha)
        y.append(val)


    print np.mean(ydata)
    print np.min(ydata)
    print np.max(ydata)

    print np.mean(xdata)
    print np.min(xdata)
    print np.max(xdata)

    plt.plot(xdata, y, 'r', label='Fit')
    plt.xlabel('log( EW [$\AA$] )')
    plt.ylabel('log ( n(EW) )')
    plt.yscale('log')
    plt.xscale('log')
    a,b = plt.ylim()
#    print a, b
#    plt.ylim([1e-5, 1e1])


    plt.xlim([1e-3, 1e2])
    plt.ylim([1e-5, 1e2])
    plt.legend(frameon=False, loc='lower left', prop={'size':8})

plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.suptitle(r'{0:s}, a={1:s}, Rvir={2:.1f} kpc'.format(galID, expn, rvir))
s = '{0:s}_{1:s}_ewdist.pdf'.format(galID, expn)
plt.savefig(s, bbox_inches='tight')






