
from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.ticker as ticker
pd.options.mode.chained_assignment = None
import scipy.optimize as opt

def mkHist(ax,x,y,z,stat,xlabel,ylabel):
    numbins = 50
    #stat = 'count'
    if 'z' in ylabel:
        binrange = [[-3,3],[0,6]]
    else:
        binrange = [[-3,3],[-3,3]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,z,
                              statistic=stat, bins=numbins,
                                 range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    h = np.log10(h)
 
    mesh = ax.pcolormesh(xedges,yedges,h)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = plt.colorbar(mesh,ax=ax,use_gridspec=True,
                    format=ticker.FuncFormatter(fmt))
    return h,xedges,yedges,binnumber



kpc2km = 3.086e16 
pc2cm = 3.086e18 
mH = 1.6737e-24 
mSun = 1.989e33 
s2yr = 3.154e7 

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'
a = 0.490

expns = np.arange(0.200,0.500,0.01)

columns = ['a m25 rho25Rvir rho25kpc rho25cokpc '
           'm50 rho50Rvir rho50kpc rho50cokpc '
           'm75 rho75Rvir rho75kpc rho75cokpc '
           'm90 rho90Rvir rho90kpc rho90cokpc'][0].split()

lowestN,highestN = -5.5,-2.5
numNbins = 6
denseBins = np.linspace(lowestN,highestN,numNbins+1)
denseLabels = []
for i in range(numNbins):
    s = 'nLim=({0:.1f},{1:.1f})'.format(denseBins[i],denseBins[i+1])
    denseLabels.append(s)

results = np.zeros((numNbins,len(expns),len(columns)))

fullCols = 'expn fullMass meanDense medianDense'.split()
fullResults = np.zeros((len(expns),len(fullCols)))
for i,a in enumerate(expns):

    print(a)

    rname = dataloc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = dataloc + filename.format(a)
    df = pd.read_hdf(fname,'data')
    df['mass'] = df['density']*mH*(df['cell_size']*pc2cm)**3 / mSun

    loT,hiT = 3.5,4.5
    tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    spaceInds = (df['theta']<80) & (df['r']>0.5*rvir)
    fullDenseInds = (df['density']>10**lowestN) & (df['density']<10**highestN)
    fullFilament = df[tempInds & spaceInds & fullDenseInds]
    totalMass = fullFilament['mass'].sum()

    fullResults[i][0] = a
    fullResults[i][1] = np.log10(totalMass)
    fullResults[i][2] = np.log10(fullFilament['density'].mean())
    fullResults[i][3] = np.log10(fullFilament['density'].median())

    for j in range(numNbins):
        
        loN = denseBins[j]
        hiN = denseBins[j+1]
    
        denseInds = (df['density']>10**loN) & (df['density']<10**hiN)

        fil = df[tempInds & denseInds & spaceInds]
        fil['rho'] = np.sqrt(fil['xRot']**2 + fil['yRot']**2)/rvir

        rhoMin = 0
        rhoMax = 3
        rhoBins = np.linspace(rhoMin,rhoMax,5000)

        mIn = []
        for k in range(len(rhoBins)):
            inside = fil['rho']<=rhoBins[k]
            massInside = fil['mass'][inside].sum()
            mIn.append(massInside)

        results[j,i,0] = a

        fractions = [0.25,0.50,0.75,0.90]
        ind = 0
        for fraction in fractions:
            m90 = mIn[-1]*fraction
            percentile = np.digitize(np.array(m90),mIn)
            containingRadius = rhoBins[percentile]

            volume = fil['r'].max()*np.pi*(containingRadius*rvir)**2

            results[j,i,ind+1] = m90/totalMass
            results[j,i,ind+2] = containingRadius
            results[j,i,ind+3] = containingRadius*rvir
            results[j,i,ind+4] = mIn[-1]/volume
            #results[j,i,ind+4] = (containingRadius*rvir)/a
            ind += 4

store = pd.HDFStore('massContained.h5',mode='w')
for i in range(numNbins):
    res = results[i]
    df = pd.DataFrame(res,columns=columns)
    store[denseLabels[i]] = df
store.close()

df = pd.DataFrame(fullResults,columns=fullCols)
df.to_hdf('totalFilamentProps.h5','data',mode='w')

