from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white')
pd.options.mode.chained_assignment = None

loc = '/Users/jacob/research/velas/vela2b/'
loc = '/mnt/cluster/abs/cgm/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i90/{2:s}/'
filename = 'vela2b-{0:d}.{1:.3f}.{2:s}.i90.abs_cells.h5'
linesname = 'lines.dat'

galNums = range(21,30)
expns = np.arange(0.200,0.550,0.01)
ions = 'HI MgII CIV OVI'.split()

fields = 'loc.spread loc.sep length numCells'.split()
props = 'vlos t n Z'.split()
for prop in props:
    fields.append('{0:s}.mean'.format(prop))
    fields.append('{0:s}.std'.format(prop))


header = [ions, fields]
header = pd.MultiIndex.from_product(header,names=['ion','field'])
index = [galNums,expnLabels]
index = pd.MultiIndex.from_product(index,names=['galNum','expn'])
results = pd.DataFrame(columns=header,index=index)


for galNum in galNums:
    for i,(a,aLabel) in enumerate(zip(expns,expnLabels)):
        for ion in ions:
            print(galNum,aLabel,ion)
            floc = loc+subloc.format(galNum,a,ion)
            fname = floc+filename.format(galNum,a,ion)
            lname = floc+'lines.dat'
            try: 
                df = pd.read_hdf(fname,'data')
                
                lheader = 'xen yen zen xex yex zex'.split()
                iheader = 'los impact phi inc'.split()
                lines = pd.read_csv(lname,sep='\s+',
                                    skiprows=2,names=lheader)
                losinfo = pd.read_csv(lname.replace('dat','info'),
                                      sep='\s+',skiprows=2,names=iheader)
                losList = df['LOS'].unique()
                spread,sep,lens,nCells = [],[],[],[]
                temp,dense,metal,vlos = [],[],[],[]
                tempStd,denseStd,metalStd,vlosStd = [],[],[],[]
                print(len(losList))
                for j,los in enumerate(losList):
                    cells = df[df['LOS']==los]
                    points = lines.iloc[int(los-1)]
                    xen,yen,zen = points[0],points[1],points[2]
                    xex,yex,zex = points[3],points[4],points[5]
                    length = np.sqrt( (xex-xen)**2 + (yex-yen)**2 + (zex-zen)**2)
                    dx = (xex-xen)/length
                    dy = (yex-yen)/length
                    dz = (zex-zen)/length
                    
                    cells['dist'] = np.sqrt((cells['x']-xen)**2 + 
                                            (cells['y']-yen)**2 +
                                            (cells['z']-zen)**2) / length
                    cells['vlos'] = (cells['vx']*dx + cells['vy']*dy + 
                                        cells['vz']*dz)
                    
                    minmax = cells['dist'].max()-cells['dist'].min()
                    #seperate = cells['dist'][1:] - cells['dist'][:-1]
                    cells['t'] = 10**cells['temperature']
                    cells['n'] = 10**cells['nH']
                    
                    lens.append(length)
                    spread.append(minmax)
                    temp.append(cells['t'].mean())
                    dense.append(cells['n'].mean())
                    metal.append(cells['alpha_Zmet'].mean())
                    vlos.append(cells['vlos'].mean())
                    nCells.append(len(cells))
                    
                    if len(cells)>1:
                        sep.append(cells['dist'].sort_values().diff().max())
                        tempStd.append(cells['t'].std())
                        denseStd.append(cells['n'].std())
                        metalStd.append(cells['alpha_Zmet'].std())    
                        vlosStd.append(cells['vlos'].std())
                                    
                results[ion,'loc.spread'].loc[galNum,aLabel] = np.mean(spread)
                results[ion,'loc.sep'].loc[galNum,aLabel] = np.mean(sep)
                results[ion,'vlos.mean'].loc[galNum,aLabel] = np.mean(vlos)
                results[ion,'vlos.std'].loc[galNum,aLabel] = np.mean(vlosStd)
                results[ion,'t.mean'].loc[galNum,aLabel] = np.mean(temp)
                results[ion,'t.std'].loc[galNum,aLabel] = np.mean(tempStd)
                results[ion,'n.mean'].loc[galNum,aLabel] = np.mean(dense)
                results[ion,'n.std'].loc[galNum,aLabel] = np.mean(denseStd)
                results[ion,'Z.mean'].loc[galNum,aLabel] = np.mean(metal)
                results[ion,'Z.std'].loc[galNum,aLabel] = np.mean(metalStd)
                results[ion,'length'].loc[galNum,aLabel] = np.mean(lens)
                results[ion,'numCells'].loc[galNum,aLabel] = np.mean(nCells) 

            except IOError:
                continue
                                    
print('Done')

outfile = 'vpAssumptions_vela.h5'
results.to_hdf(outfile,'data')




