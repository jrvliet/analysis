#!/usr/bin/python

# Cacluates the slope of the EW vs. D results from
# master_bulk

import numpy as np

galID_list = ['D9o2', 'D9q', 'D9m4a']
labels = ['dwSN', 'dwALL\_1', 'dwALL\_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']

ewcut = 0

for ion in ion_list:

    for galID in galID_list:

        outfile = galID+'.'+ion+'.cut_0mA_ewImp.out'

        data = np.loadtxt(outfile, skiprows=1)

        imp = data[:,0] 
        ew = data[:,1]
        err = data[:,2]

        print ew

        ew = np.log10(ew)
        err = np.log10(err)
        w = [1/(i*i) for i in err]

        p, r = np.polyfit(imp, ew, 1)

        print ion, galID, p, r
