from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None

loc = '/mnt/cluster/abs/cgm/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i90/{2:s}/'
filename = 'vela2b-{0:d}.{1:.3f}.{2:s}.i90.abs_cells.h5'
linesname = 'lines.dat'
ions = 'HI MgII CIV OVI'.split()

galNums = range(21,30)
expns = np.arange(0.200,0.500,0.01)

header = 'expn redshift'.split()
for ion in ions:
    header.append('{0:s}.mean'.format(ion))
    header.append('{0:s}.std'.format(ion))


for galNum in galNums:
    print(galNum)
    spread = np.zeros((len(expns),len(header)))
    spread = pd.DataFrame(spread,columns=header)
    for i,a in enumerate(expns):
        spread['expn'].iloc[i] = a
        spread['redshift'].iloc[i] = 1./a - 1
        for ion in ions:
            l = loc+subloc.format(galNum,a,ion)
            fname = l+filename.format(galNum,a,ion)
            lname = l+linesname.format(ion)
            try:
                df = pd.read_hdf(fname,'data')
                lheader = 'xen yen zen xex yex zex'.split()
                iheader = 'los impact phi inc'.split()
                lines = pd.read_csv(lname,sep='\s+',
                                    skiprows=2,names=lheader)
                losinfo = pd.read_csv(lname.replace('dat','info'),
                                      sep='\s+',skiprows=2,names=iheader)

                losList = df['LOS'].unique() 
                s = []
                for j,los in enumerate(losList):
                    cells = df[df['LOS']==los]
                    points = lines.iloc[int(los-1)]
                    impact = losinfo['impact'].iloc[int(los-1)]
                    xen,yen,zen = points[0],points[1],points[2]
                    xex,yex,zex = points[3],points[4],points[5]
                    length = np.sqrt( (xex-xen)**2 + (yex-yen)**2 + (zex-zen)**2)
                    cells['dist'] = np.sqrt((cells['x']-xen)**2 + 
                                            (cells['y']-yen)**2 +
                                            (cells['z']-zen)**2) / length
                    minmax = cells['dist'].max()-cells['dist'].min()
                    #print(los,impact,minmax)
                    s.append(minmax)
                spread[ion+'.mean'].iloc[i] = np.mean(s)
                spread[ion+'.std'].iloc[i] = np.std(s)
            except IOError:
                print('File not found: {0:s}'.format(fname))
                continue
    outfile = 'vela2b-{0:d}_absorberLocSpread.h5'.format(galNum)
    spread.to_hdf(outfile,'data')
