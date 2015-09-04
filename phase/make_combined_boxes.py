#!/usr/bin/python

import numpy as np


galID_list = ['D9o2', 'D9q', 'D9m4a']
galID_list = ['D9o2']
feedback_list = ['dwSN', 'dwALL_1', 'dwALL_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]

gasfile_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/gasfiles/'


for i in range(0, len(galID_list)):

    galID = galID_list[i]
    expn = expn_list[i]

    for ion in ion_list:

        combined_f = galID+'_'+ion+'_GZ_combined.txt'
        comf = open(gasfile_loc+combined_f, 'w')

        for a in expn:
            
            gasfile = galID+'_GZa'+a+'.'+ion+'.txt'
            gasf = open(gasfile_loc+gasfile)
            gasf.readline()
            header = gasf.readline()

            if expn.index(a)==0:
                comf.write(header)
                
            for line in gasf:
                comf.write(line)

            gasf.close()

        comf.close()
        
