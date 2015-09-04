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

if len(sys.argv)!=4:
    print '\nUsage:'
    print 'genericEWvS <galID> <expn> <Rvir>'
    print ''
    sys.exit()
galID = sys.argv[1]
expn = sys.argv[2]
rvir = float(sys.argv[3])


ion_list = ['HI', 'MgII', 'CIV', 'OVI']
#galID_list = ['D9o', 'D9o2']
#a_list = ['1.001', '1.002']
#c = ['blue', 'red']

Dmin =  []
Dmax = []
for i in range(0,15):
    Dmin.append(i*0.1)
    Dmax.append((i+1)*0.1)
    
for ion in ion_list:

#    filename = './'+ion+'/'+galID+'_GZa'+expn+'.'+ion+'.txt'
    filename = './'+ion+'/'+galID+'.'+ion+'.a'+expn+'.ALL.sysabs'
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


    plt.errorbar(imp_mean, ew_mean, yerr=ew_err, xerr=[imp_err_neg, imp_err_pos], fmt='s', linewidth=2)
    plt.title(ion)
    plt.yscale('log')
    plt.xlim([0,1.5])
    plt.xlabel('Impact Parameter [Rvir]')
    plt.ylabel('log(EW)')
    plt.savefig(galID+'_'+ion+'_EWvsD.pdf')
    plt.clf()
    plt.cla()
