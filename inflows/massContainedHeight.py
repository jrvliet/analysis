
'''
Determines the width of the filament by looking at the radius that contains 90%
of the mass. Cacluates how this changes with distance from the galaxy
'''


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

store = pd.HDFStore('massContainedHeight.h5',mode='w')

expns = np.arange(0.200,0.550,0.01)

header = 'loz hiz meanz rho90rvir rho90pkpc massFrac mass numcells'.split()

# Phase selection
loT,hiT = 3.5,4.5
lowestN,highestN = -5.5,-2.5
numNbins = 3
denseBins = np.linspace(lowestN,highestN,numNbins+1)
denseLabels = 'low mid high'.split()

minHeight, maxHeight = 0,6
numHeightBins = 12
heightBins = np.linspace(minHeight,maxHeight,numHeightBins+1)

for i,a in enumerate(expns):

    print(a)

    rname = dataloc+rotname.format(a)
    with open(rname) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    fname = dataloc + filename.format(a)
    df = pd.read_hdf(fname,'data')
    df['mass'] = df['density']*mH*(df['cell_size']*pc2cm)**3 / mSun

    tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    spaceInds = (df['theta']<80) & (df['r']>0.5*rvir)
    fullDenseInds = (df['density']>10**lowestN) & (df['density']<10**highestN)
    fullFilament = df[tempInds & spaceInds & fullDenseInds]
    totalMass = fullFilament['mass'].sum()

    for j in range(numNbins):
        
        loN = denseBins[j]
        hiN = denseBins[j+1]
        print('\tloN,hiN = {0:.1f},{1:.1f}'.format(loN,hiN))
    
        denseInds = (df['density']>10**loN) & (df['density']<10**hiN)

        results = np.zeros((numHeightBins,len(header)))
    
        for k in range(numHeightBins):

            # Select out the height slice
            loz = heightBins[k]*rvir
            hiz = heightBins[k+1]*rvir
            heightInds = (df['zRot']>loz) & (df['zRot']>hiz)
        

            # Select out filament material
            fil = df[tempInds & denseInds & spaceInds & heightInds]
            fil['rho'] = np.sqrt(fil['xRot']**2 + fil['yRot']**2)/rvir

            print('\t\tloz,hiz = {0:.1f},{1:.1f}\tnumcells = {2:d}'.format(
                    loz,hiz,len(fil)))

            if len(fil)>0:

                meanz = fil['zRot'].mean()/rvir
                rhoMin = 0
                rhoMax = 3
                rhoBins = np.linspace(rhoMin,rhoMax,5000)

                # Move through bins of rho
                # Determine mass within each distance
                mIn = []
                for l in range(len(rhoBins)):
                    insideInds = fil['rho']<=rhoBins[l]
                    massInside = fil['mass'][insideInds].sum()
                    mIn.append(massInside)
                
                fig,ax = plt.subplots(1,1,figsize=(5,5))
                ax.plot(rhoBins,mIn)
                ax.set_xlabel('Rho')
                ax.set_ylabel('mIn')
                fig.savefig('mIn.png')
                plt.close(fig)

                totalMass = mIn[-1]
                # Determine the distance that contains 90% of the mass
                fraction = 0.90
                m90 = mIn[-1]*fraction
                percentile = np.digitize(np.array(m90),mIn)
                containingRadius = rhoBins[percentile]

                volume = fil['r'].max()*np.pi*(containingRadius*rvir)**2

            else:
                meanz = np.nan
                containingRadius = np.nan
                m90 = np.nan
                totalMass = np.nan
    
            results[k,0] = loz/rvir
            results[k,1] = hiz/rvir
            results[k,2] = meanz
            results[k,3] = containingRadius
            results[k,4] = containingRadius*rvir
            results[k,5] = m90/totalMass
            results[k,6] = np.log10(totalMass)
            results[k,7] = len(fil)
                
        denseLabel = denseLabels[j]
        expnLabel = 'a{0:d}'.format(int(a*1000))
        path = '{0:s}/{1:s}'.format(expnLabel,denseLabel)
        resultsdf = pd.DataFrame(results,columns=header)

        store[path] = resultsdf

store.close()
