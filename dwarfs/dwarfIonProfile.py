'''
Plots the radial profile of each species and it's ionization fraction in aggregate
to compare the difference feedback regiemes
'''

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

loc = '/lustre/projects/p089_swin/simulations/dwarfs/'
filename = '{0:s}_outputs/transfer/{0:s}_GZa*.{1:s}.txt'

loc = '/home/jacob/research/dwarfs/gasfiles/'
filename = '{0:s}_GZa*.{1:s}.txt'

galIDs = 'D9o2 D9q D9m4a'.split()
galNames = 'dwSN dwALL1 dwALL8'.split()
ions = 'HI MgII CIV OVI'.split()
elements = 'H Mg C O'.split()
header = ['size x y z vx vy vz nH temperature SNII SNIa '
          'nAtom fIon nIon alphaSol alphaZmet id '
          'tph trec tcoll tcool']
header = header[0].split()
numbins = 20

fields = 'nH temperature alphaSol nAtom fIon nIon'.split()
labels = 'nH T Z $n_{0:s}$ $f_{0:s}$ $n_{0:s}$'.split()
for ion,e in zip(ions,elements):

    fig,axes = plt.subplots(3,2,figsize=(10,15))

    for galID,galName in zip(galIDs,galNames):

        fpattern = loc+filename.format(galID,ion)
        files = glob.glob(fpattern)

        print(fpattern)
        gas = pd.read_csv(files[0],sep='\s+',names=header,skiprows=2)

        for fname in files[1:]:
            print(fname)
            df = pd.read_csv(fname,sep='\s+')
            gas = gas.append(df,ignore_index=True)

        rvir = gas['x'].max()/2.
        gas['r'] = np.sqrt(gas['x']**2 + gas['y']**2 + gas['z']**2)/rvir

        # Bin the data
        bins = np.linspace(gas['r'].min(),gas['r'].max(),numbins+1)
        gas['binned'] = pd.cut(gas['r'],bins,labels=range(numbins))
        values = gas.groupby(gas['binned']).mean()
        errs = gas.groupby(gas['binned']).std()

        for i,(ax,field,label) in enumerate(zip(axes.flatten(),fields,labels)):
        
            ax.errorbar(values['r'],values[field],yerr=errs[field],label=galName)
            ax.set_xlabel('Distance [Rvir]')
            if i<3:
                l = label
            elif i==3:
                l = label.format(e)
            else:
                l = label.format(ion)
            ax.set_ylabel('log( {0:s} )'.format(l))

    axes[0,0].legend(loc='best')
    fig.tight_layout()
    s = 'dwarfProfile_{0:s}.pdf'.format(ion)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)
        

            
            





