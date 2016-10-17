
'''
Creates narrow density cuts through the inflow
Fit a rayleigh distribution to the distribution of
distance from each cell to the line of best fit.
Also calculates the sigma of the data and FWHM of 
fit. Plots vs time.
'''


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import scipy.stats as st
import scipy.linalg as sl
import scipy.interpolate as sint

pd.options.mode.chained_assignment = None

def fit_line(locM):
    
    covar = locM.cov()
    (u,s,v) = sl.svd(covar)
    u0 = u[0,0]
    u1 = u[1,0]
    u2 = u[2,0]
    
    return u0,u1,u2

def get_fwhm(param,xmin,xmax):

    xmin = 0
    x = np.linspace(xmin,xmax,10000)
    y = st.rayleigh.pdf(x,loc=param[0],scale=param[1])

    halfmax = max(y)/2.

    spline = sint.UnivariateSpline(x,y-halfmax, s=0)
    r = spline.roots()
    
    return r



dataloc = '/home/jacob/research/velas/vela2b/vela27/'
dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotmat = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expns = np.arange(0.200,0.500,0.01)

# Cloud limits
loT, hiT = 10**3.5, 10**4.5

# Density bins
nBins = np.linspace(-5.5,-2.5,7)

header = ['a','redshift','loN','hiN','numCells','stdDev','rayleighLoc',
            'rayleighScale','rayleighfwhm','rayleighStd','speedStd',
            'valongStd','vperpStd','xRotStd','yRotStd','zRotStd','zRotRange',
            'snIIStd','snIaStd','snIImean','snIamean','rMean','rStd','nHmean',
            'nHStd','tMean','tStd','rMeanMod','speedMean','valongMean','vperpMean']
fit = np.zeros(len(header))

for a in expns:
    
    print(a)

    # Read in data
    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname, 'data')

    # Select regions
    tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
    spacInd = (df['x']<0) & (df['z']>0)

    # Get Rvir
    with open(dataloc+rotmat.format(a)) as f:
        f.readline()
        rvir = float(f.readline().split()[3])

    # Loop over density cuts 
    nCombs = it.combinations(nBins,2)
    for combo in nCombs:
        if combo[1]>combo[0]:

            loN = combo[0]
            hiN = combo[1]
            densInd = (df['density']>10**loN) & (df['density']<10**hiN)
            thisfit = np.zeros(len(header))
            thisfit[0] = a
            thisfit[1] = 1./a-1
            thisfit[2] = loN
            thisfit[3] = hiN
        
            cloud = df[densInd & spacInd & tempInd]
            thisfit[4] = len(cloud)

            cloud['speed'] = np.sqrt(cloud['vx']**2 + cloud['vy']**2 + cloud['vz']**2)

            if len(cloud)>100:
    
                cloud['r'] = np.sqrt(cloud['x']**2 + cloud['y']**2 + cloud['z']**2)
                # Fit line
                cloudLoc = cloud[['x','y','z']]
                locM = cloudLoc - cloudLoc.mean()
                u0,u1,u2 = fit_line(locM)
            
                # Get distance from each cell to this line
                locM['a0'] = u2*locM['y'] - u1*locM['z']
                locM['a1'] = u0*locM['z'] - u2*locM['x']
                locM['a2'] = u1*locM['x'] - u0*locM['y']
                locM['dist'] = np.sqrt(locM['a0']**2 + locM['a1']**2 + locM['a2']**2)
                
                # Fit Rayleigh to distribution
                param = st.rayleigh.fit(locM['dist'])
                roots = get_fwhm(param,locM['dist'].min(),locM['dist'].max())
                if len(roots)==2:
                    fwhm = roots[1]-roots[0] 
                else:
                    fwhm = np.NAN
                    
                # Get kinematics
                adotb = cloud['vx']*u0 + cloud['vy']*u1 + cloud['vz']*u2
                bdotb = u0**2 + u1**2 + u2**2
                factor = adotb/bdotb
                cloud['along0'] = factor*u0
                cloud['along1'] = factor*u1
                cloud['along2'] = factor*u2
                cloud['along'] = np.sqrt(cloud['along0']**2 +
                                cloud['along1']**2 + cloud['along2']**2)
                cloud['perp0'] = cloud['vx'] - cloud['along0']
                cloud['perp1'] = cloud['vy'] - cloud['along1']
                cloud['perp2'] = cloud['vz'] - cloud['along2']
                cloud['perp'] = np.sqrt(cloud['perp0']**2 +
                                cloud['perp1']**2 + cloud['perp2']**2)
                
                # Fill output array
                thisfit[5] = locM['dist'].std()
                thisfit[6] = param[0]
                thisfit[7] = param[1]
                thisfit[8] = fwhm
                thisfit[9] = st.rayleigh.std(loc=param[0],scale=param[1])
                thisfit[10] = cloud['speed'].std()
                thisfit[11] = cloud['along'].std()
                thisfit[12] = cloud['perp'].std()
                thisfit[13] = cloud['xRot'].std()
                thisfit[14] = cloud['yRot'].std()
                thisfit[15] = cloud['zRot'].std()
                thisfit[16] = cloud['zRot'].max() - cloud['zRot'].min()
                thisfit[17] = cloud['SNII'].std()
                thisfit[18] = cloud['SNIa'].std()
                thisfit[19] = cloud['SNII'].mean()
                thisfit[20] = cloud['SNIa'].mean()
                thisfit[21] = cloud['r'].mean()
                thisfit[22] = cloud['r'].std()
                thisfit[23] = cloud['density'].mean()
                thisfit[24] = cloud['density'].std()
                thisfit[25] = cloud['temperature'].mean()
                thisfit[26] = cloud['temperature'].std()
                thisfit[27] = cloud['r'].mean()/rvir
                thisfit[28] = cloud['speed'].mean()
                thisfit[29] = cloud['along'].mean()
                thisfit[30] = cloud['perp'].mean()
            
            else:
                thisfit[5:] = np.NAN
                thisfit[10] = cloud['speed'].std()
                
            fit = np.vstack((fit,thisfit))
    

fit = np.delete(fit,(0),axis=0)
df = pd.DataFrame(fit,columns=header)
outfile = 'projected_distance_distribution.h5'
df.to_hdf(outfile,'data',mode='w')

            

















