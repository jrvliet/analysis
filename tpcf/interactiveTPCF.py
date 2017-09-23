
'''
Code that allows for custom TPCF samples.
Mixes samples from muptiple galaxies, inclinations, azimuthal angles
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import itertools as it
import joblib as jl
import numb as nb
import tempfile
import sys
import decimal
import os


class tpcfProps(object):
    '''
    Class to define settings for TPCF
    '''
    
    def __init__ (self):
        self.ewLo = 0.
        self.ewHi = 10.
        self.dLo = 0.
        self.dHi = 200.
        self.azLo = 0.
        self.azHi = 90.
        self.binSize = 10.
        self.bootNum = 1000


class tpcfRun(object):
    '''
    Class that holds the settings for this interactive run
    '''
    
    def __init__ (object):
        self.azLo = 0.
        self.azHi = 90.
        self.expn = 0.490
        self.ions = 'MgII CIV OVI'.split()
        self.iLo = 0.
        self.iHi = 90.
        self.dLo = 0.
        self.dHi = 200.
        self.loc = '/mnt/cluster/abs/cgm/vela2b/'
        
def read_input():
    '''
    Reads in the input file and fills out run object
    '''
    
    fname = 'tpcf.config'
    run = tpcfRun()
    with open(fname,'r') as f:
        
        run.azLo = float(f.readline().split()[0])
        run.azHi = float(f.readline().split()[0])

        run.dLo = float(f.readline().split()[0])
        run.dHi = float(f.readline().split()[0])

        run.iLo = float(f.readline().split()[0])
        run.iHi = float(f.readline().split()[0])
        run.expn = float(f.readline().split()[0])

        # Read in ions
        f.readline()
        ions = []
        for line in f:
            ions.append(line.split()[0])
        run.ions = ions
            
    return run
        

def select_los(run):
    '''
    Selects lines of sight that fit the limits contained in run
    '''

    linesHeader = 'los impact phi incline az'.split()
    galNums = range(21,30)
    los = pd.dataframe(columns=galNums)
    for galNum in galNums:
        
        # Check if the expansion parameter exists
        dirname = run.loc+'vela{0:d}/a{1:.3f}'.format(galNum,run.expn)
        if os.path.isdir(dirname):
            
            # Get list of inclinations in this directory
            iDirs = [name for name in os.listdir('.') if 
                        os.path.isdir(os.path.join('.',name)) and
                        i[0]=='i' and 
                        float(i.split('i')[1])>=run.iLo and
                        float(i.split('i')[1])<=run.iLo]
            for iDir in iDirs:
                linesFile = '{0:s}/{1:s}/lines.info'.format(dirname,iDir)
                lines =  pd.read_csv(linesFile,names=linesHeader,
                                     sep='\s+')
                targets = ((lines['az']>=run.azLo) & 
                           (lines['az']<=run.azHi) &
                           (lines['impact']>=run.dLo) & 
                           (lines['impact']<=run.dHi))
                target = lines['los'][targets]
                
                los.set_value(iDir,galNum,target)
            
        
    
        # Read in the lines file
    
        # Select los that fit the requirements

        # Add to array
        
    # Return array

    













