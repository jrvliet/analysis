
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
import numba as nb
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
    
    def __init__ (self):
        self.azLo = 0.
        self.azHi = 90.
        self.expn = '0.490'
        self.ions = 'MgII CIV OVI'.split()
        self.iLo = 0.
        self.iHi = 90.
        self.dLo = 0.
        self.dHi = 200.
        self.loc = '/mnt/cluster/abs/cgm/vela2b/'

    def print_run(self):
        print('\nTPCF Run Properties: ')
        print('\tLocation = {0:s}'.format(self.loc))
        print('\tExpansion Parameter = {0:s}'.format(self.expn))
        print('\tIons = {0:s}'.format(', '.join(self.ions)))
        print('\tAzimuthal Limits = {0:f} - {1:f}'.format(self.azLo,self.azHi))
        print('\tInclination Limits = {0:f} - {1:f}'.format(self.iLo,self.iHi))
        print('\tImpact Limits = {0:f} - {1:f}'.format(self.dLo,self.dHi))
        print()
        
def read_input():
    '''
    Reads in the input file and fills out run object
    '''
    
    fname = 'tpcf.config'
    run = tpcfRun()
    with open(fname,'r') as f:
        
        run.expn = f.readline().split()[0]

        run.azLo = float(f.readline().split()[0])
        run.azHi = float(f.readline().split()[0])

        run.dLo = float(f.readline().split()[0])
        run.dHi = float(f.readline().split()[0])

        run.iLo = float(f.readline().split()[0])
        run.iHi = float(f.readline().split()[0])

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
    iDirs,iDirsList = find_inclinations(run,galNums)
        
    los = pd.DataFrame(columns=pd.MultiIndex.from_tuples(iDirs),
                        index=range(1000))
    
    # Loop through to read in lines.info files
    for gal,inc in iDirsList:
        dirname = '{0:s}/vela{1:d}/{2:s}/{3:s}/'.format(run.loc,
                                                gal,run.expn,inc)
        linesFile = '{0:s}/lines.info'.format(dirname)
        lines =  pd.read_csv(linesFile,names=linesHeader,
                             sep='\s+')
        targets = ((lines['az']>=run.azLo) & 
                   (lines['az']<=run.azHi) &
                   (lines['impact']>=run.dLo) & 
                   (lines['impact']<=run.dHi))
        target = lines['los'][targets]
        
        los[gal,inc] = pd.Series(target)
        los.set_value(iDir,galNum,target)
        
    return los
   
def build_sample(run,los):
    '''
    Using LOS contained in los, select out the velocity differences 
    from TPCF files
    '''
    
    # Loop through columns in los
    allVels = [pd.DataFrame for ion in ions]
    for galNum,inc in los.columns:
    
        dirname = '{0:s}/vela{1:d}/{2:s}/{3:s}/tpcf/'.format(run.loc,
                    galNum,run.expn,inc)
        for i,ion in eunumerate(run.ions):
            filename = 'vela2b-{0:d}_{1:s}_{2:s}_{3:s}_velDiff.csv'.format(
                    galNum,run.expn,inc,ion)
            fname = dirname+filename
            vd = pd.read_csv(fname)
            velDiffColumns = set(vd.columns.values)
            
            # Select out the lines that are in LOS
            snap = los[galNum,inc]
            a = set(['{0:.1f}'.format(v) for v in snap.values])
            snapVels = vd[list(a & velDiffColumns)]
            allVels[i] = pd.concat([allVels[i],snapVels],axis=1)
    
    # Reset allVels column names
    for df in allVels:
        df.columns = range(df.shape[1])

    # Convert to memmap objects
    allVelsPaths = []
    allVelsShapes = []
    maxVel = 0
    for df,ion in zip(allVels,ions):
        path = tempfile.mkdtemp()
        velMemPath = os.path.join(velDiffPath,
                    'vellDiff_{0:s}.mmap'.format(ion))
        velDiffMem = np.memmap(velMemPath,dtype='float',
                        shape=df.shape,mode='w+')
        velDiffMem[:] = df.values[:]
        allVelsPaths.append(velMemPath)
        dfMax = np.nanmax(df.values)
        if dfMax>maxVel:
            maxVel = dfMax

    return allVelsPath,allVelsShapes,maxVel
    

def sample_bins(run,maxVel,tpcfProp):
    ''' 
    Generates the velocity bins and labels to making the tpcf
    '''

    nbins = int(maxVel/tpcfProp.binSize)
    endPoint = tpcf.binSize*(nbins+1)
    bins = np.arange(0,endPoint,tpcf.binSize)

    labels = [(bins[i]+bins[i+1])/2. for i in range(nbins)]
    lastLabel = labels[-1] + (labels[1]-labels[0])
    labels.append(lastLabel)
    return bins,labels
    


def sample_tpcf(run,samplePath,sampleShape,bins,labels,bootstrap=0):
    '''
    Constructs the TPCF from the sample
    '''
    
    sample = np.memmap(samplePath,dtype='float',mode='r',
                        shape=sampleShape)
    if bootstrap!=0:
        sample = sample[:,np.random.random.choice(sample.shape[1],
                        sample.shape[1],replace=True)]
    
    flat = sample.flatten()
    flat = flat[~np.isnan(flat)]
    tpcf = np.sort(np.bincount(np.digitize(flat,bins)))[::-1]
    tpcf = tpcf/tpcf.sum()
    return tpcf 
    
    
    
    

def find_inclinations(run,galNums):

    '''
    Returns a list of inclinations directories for each galaxy number
    '''

    iDirs = {}
    
    for galNum in galNums:

        # Check if the expansion parameter exists
        dirname = run.loc+'vela{0:d}/a{1:.3f}'.format(galNum,run.expn)
        inclines = []
        if os.path.isdir(dirname):
            
            # Get list of inclinations in this directory
            inclines = [name for name in os.listdir('.') if 
                        os.path.isdir(os.path.join('.',name)) and
                        names[0]=='i' and 
                        float(names.split('i')[1])>=run.iLo and
                        float(names.split('i')[1])<=run.iHi]

        iDirs[galNum] = inclines

    iDirsList = []
    for key,val in iDirs.items():
        for v in val:
            iDirsList.append([key,v])
    
    return iDirs,iDirsList
    



if __name__ == '__main__':

    run = read_input()
    run.print_run()
    tpcfProp = tpcfProps()
    tpcfProp.bootNum = 10
    
    los = select_los(run)
    allVelsPath,allVelsShapes,maxVel = build_sample(run,los):
    








