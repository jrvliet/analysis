#!/usr/bin/python


"""
Usage:
python genericEWvD <galID> <expn> <Rvir>

import numpy as np
import matplotlib.pyplot as plt
"""
import sys
import matplotlib.pyplot as plt
import numpy as np


# Read in galaxy properties from galaxy.props
f = open('galaxy.props')
galID = f.readline().split()[1]
expn = f.readline().split()[1]
redshift = f.readline().split()[1]
mvir = float(f.readline().split()[1])
rvir = float(f.readline().split()[1])
inc = int(float(f.readline().split()[1]))




ion_list = ['HI', 'MgII', 'CIV', 'OVI']
#galID_list = ['D9o', 'D9o2']
#a_list = ['1.001', '1.002']
#c = ['blue', 'red']

Dmin =  []
Dmax = []
for i in range(0,15):
    Dmin.append(i*0.1)
    Dmax.append((i+1)*0.1)
    
ewAll = []
impAll = []
yerrAll = []
xerrPosAll = []
xerrNegAll = []
for ion in ion_list:

#    filename = './'+ion+'/'+galID+'_GZa'+expn+'.'+ion+'.txt'
#    filename = './'+ion+'/'+galID+'.'+ion+'.a'+expn+'.ALL.sysabs'
    filename = './{0:s}/{1:s}.{0:s}.a{2:s}.i{3:d}.ALL.sysabs'.format(ion,galID,expn,inc)
    impRaw, ewRaw = np.loadtxt(filename, skiprows=1, usecols=(1, 5), unpack=True)
#    ew1_raw, imp1_raw = np.loadtxt(fname1, skiprows=1, usecols=(5, 22), unpack=True)
#    ew2_raw, imp2_raw = np.loadtxt(fname2, skiprows=1, usecols=(5, 22), unpack=True)

    # Remove any EW that is zero, as these come from LOS with
    # no detections
    ew = []
    imp = []
    for i in range(0,len(ewRaw)):
        if ewRaw[i]!=0.0:
            ew.append(ewRaw[i])
            imp.append(impRaw[i]/rvir)
            
    ew = np.array(ew)
    imp = np.array(imp)

#    print len(ew)
#    print ew
#    print imp
#    print np.mean(ew)


    ew_mean = []
    ew_err = []
    imp_mean = []
    imp_err_pos = []
    imp_err_neg = []
    for i in range(0,len(Dmin)):
        
        low = Dmin[i]
        high = Dmax[i]
        
        ewtotal = []
        imptotal = []
        for j in range(0,len(imp)):
            
            if imp[j]<high and imp[j]>low:
                ewtotal.append(ew[j])
                imptotal.append(imp[j])
        
        im = np.mean(imptotal)
        imp_mean.append(im)
        ew_mean.append(np.mean(ewtotal))
        ew_err.append(np.std(ewtotal))
        imp_err_pos.append(high-im)
        imp_err_neg.append(im-low)

    ewAll.append(ew_mean)
    impAll.append(imp_mean)
    yerrAll.append(ew_err)
    xerrNegAll.append(imp_err_neg)
    xerrPosAll.append(imp_err_pos)


for i, ion in enumerate(ion_list):
    subplotnum = 221+i
    plt.subplot(subplotnum)
    plt.errorbar(impAll[i], ewAll[i], yerr=yerrAll[i], 
                 xerr=[xerrNegAll[i], xerrPosAll[i]], 
                 fmt='s', linewidth=2)

    plt.yscale('log')
    plt.xlim([0,1.5])
    plt.xlabel('Impact Parameter [Rvir]')
    plt.ylabel(r'$\log(EW_{'+ion+'})$')

plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.suptitle('{0:s}, a={1:s}, Rvir={2:.1f} kpc, i={3:d}'.format(galID, expn, rvir, inc))
name = '{0:s}_a{1:s}_i{2:d}_EWvsD.pdf'.format(galID, expn, inc)
plt.savefig(name, bbox_inches='tight')




