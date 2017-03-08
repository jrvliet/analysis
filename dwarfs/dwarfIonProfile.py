'''
Plots the radial profile of each species and it's ionization fraction in aggregate
to compare the difference feedback regiemes
'''

from __future__ import print_functions
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

loc = '/lustre/projects/p089_swin/simulations/dwarfs/'
filename = '{0:s}_outputs/transfer/{0:s}_GZa*.{1:s}.txt'
galIDs = 'D9o2 D9q D9m4a'.split()
galNames = 'dwSN dwALL1 dwALL8'.split()
ions = 'HI MgII CIV OVI'.split()
elements = 'H Mg C O'.split()
header = 'size x xy z vx vy vz nH temperature SNII SNIa nAtom fIon nIon alphaSol alphaZmet id'.split()

for ion,e in zip(ions,elements):

    fig,axes = plt.subplots(2,2,figsize=(10,10))

    for galID,galName in zip(galIDs,galNames):

        fpattern = loc+filename.format(galID,ion)
        files = glob.glob(fpattern)

        gas = pd.read_csv(files[0],sep='\s+',names=header,skiprows=2)

        for fname in files[1:]:
            df = pd.read_csv(fname,sep='\s+')
            gas = gas.append(df,ignore_index=True)


            
            





