#!/usr/bin/python


import numpy as np
import matplotlib.pyplot as plt


ion_list = ['HI', 'MgII', 'CIV', 'OVI']
galID_list = ['D9o', 'D9o2']
a_list = ['1.001', '1.002']
c = ['blue', 'red']
Dmin =  []
Dmax = []
for i in range(0,15):
    Dmin.append(i*0.1)
    Dmax.append((i+1)*0.1)
    
for ion in ion_list:

    
    fname1 = 'D9o.'+ion+'.a1.001.ALL.sysabs.large'
    fname2 = 'D9o2.a1.002.'+ion+'.ALL.sysabs.large'
    
    ew1_raw, imp1_raw = np.loadtxt(fname1, skiprows=1, usecols=(5, 22), unpack=True)
    ew2_raw, imp2_raw = np.loadtxt(fname2, skiprows=1, usecols=(5, 22), unpack=True)

    # Remove any EW that is zero, as these come from LOS with
    # no detections
    ew1 = []
    imp1 = []
    for i in range(0,len(ew1_raw)):
        if ew1_raw[i]!=0.0:
            ew1.append(ew1_raw[i])
            imp1.append(imp1_raw[i])
            
    ew1 = np.array(ew1)
    imp1 = np.array(imp1)

    # Remove any EW that is zero, as these come from LOS with
    # no detections
    ew2 = []
    imp2 = []
    for i in range(0,len(ew2_raw)):
        if ew2_raw[i]!=0.0:
            ew2.append(ew2_raw[i])
            imp2.append(imp2_raw[i])
            
    ew2 = np.array(ew2)
    imp2 = np.array(imp2)

    ewbig = [ew1, ew2]
    impbig = [imp1, imp2]

    for k in range(0,len(ewbig)):
        galID = galID_list[k]
        ew = ewbig[k]
        imp = impbig[k]
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
    

        plt.errorbar(imp_mean, ew_mean, yerr=ew_err, xerr=[imp_err_neg, imp_err_pos], color=c[k], fmt='s', linewidth=2, label=galID)


    plt.legend(frameon=False)
    plt.savefig('ew_d_D9_comp'+ion+'.pdf')
    plt.clf()
    plt.cla()
