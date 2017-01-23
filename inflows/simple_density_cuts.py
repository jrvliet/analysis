
'''
Creates narrow density cuts through the inflow
Fit a rayleigh distribution to the distribution of
distance from each cell to the line of best fit.
Also calculates the sigma of the data and FWHM of 
fit. Plots vs time.

This version uses only simple density cuts
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
import sys

pd.options.mode.chained_assignment = None

def mfToZ(Z_cell):

    # mfToZ.py
    #   
    # Converts mass fractions to Z/Z_sun
    #   
    # Z/Z_sun = (Z_m / X_h)_cell
    #           ----------------
    #           (Z_m / X_h)_sun
    #   
    # where X_h + Y_he + Z_m = 1
    # Since the simulation does not track Helium, need to assume 
    # a value for r = Y/X
    

    # Solar values
    # Taken from Chris' Notes 
    X_sun = 0.70683
    Y_sun = 0.27431
    Z_sun = 0.0188
    
    r = 0.3347

    # Loop through cells
    X_cell = (1 - Z_cell) / (1 + r)
    Z = (Z_cell / X_cell) / (Z_sun / X_sun)

    return Z

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


# Conversion factors
u2g = 1.661e-24
g2M = 1.989e33
pc2cm = 3.086e18
boltz = 1.3807e-16
cm2km = 1e5

dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotmat = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expns = np.arange(0.200,0.550,0.01)

# Cloud limits
coolTemps = [10**3.5, 10**4.5]
warmTemps = [10**4.5, 10**5.5]
hotTemps  = [10**5.5, 10**6.5]
temps = [coolTemps,warmTemps,hotTemps]
tempLabels = ['cool','warm','hot']

# Density bins
lowestN, highestN = -5.5,-2.5
numNBins = 4
nBins = np.linspace(lowestN,highestN,numNBins)

header = 'a redshift loN hiN numCellsFrac'.split()
fields = ('speed valong vperp '
          'xRot yRot zRot r rMod rRot thetaRot phiRot '
          'SNII SNIa  density temperature mass pressure '
          'vr vzRot vrhoRot vthetaRot vphiRot thermalV vrFrac '
          'vrFracAbs metallicity'.split())
fitFields = 'stdDev rayleighLoc rayleighScale rayleighfwhm rayleighStd'.split()

header = header + fitFields
for field in fields:
    header.append(field+'Mean')
    header.append(field+'Std')
    header.append(field+'Ratio')
    header.append(field+'MeanMW')
    header.append(field+'Median')

with open('denseCutHeaders.txt','w') as f:
    for i,h in enumerate(header):
        f.write('{0:d}\t{1:s}\n'.format(i,h))

excludeCGM = False

for i in range(len(temps)):

    print(tempLabels[i])
    loT, hiT = temps[i][0], temps[i][1]
    fit = np.zeros(len(header))

    for a in expns:
        
        print('\t{0:.3f}'.format(a))

        # Get Rvir
        with open(dataloc+rotmat.format(a)) as f:
            f.readline()
            rvir = float(f.readline().split()[3])

        # Read in data
        fname = dataloc+filename.format(a)
        df = pd.read_hdf(fname, 'data')

        df['r'] = np.sqrt(df['x']**2 + df['y']**2 + df['z']**2)
        df['metallicity'] = mfToZ(df['SNII']+df['SNIa'])

        # Select regions
        tempInd = (df['temperature']>loT) & (df['temperature']<hiT)
        #spacInd = (df['x']<0) & (df['z']>0)
        fullDenseInd = (df['density']>10**lowestN) & (df['density']<10**highestN)
        if excludeCGM:
            spacInd = (df['x']<0) & (df['z']>0) & (df['r']>0.5*rvir)
        else:
            spacInd = (df['x']<0) & (df['z']>0)

        fullCloudSize = float((tempInd & spacInd & fullDenseInd).sum())

        # Loop over density cuts 
        nCombs = it.combinations(nBins,2)
        for j in range(len(nBins)-1):
            loN,hiN = nBins[j],nBins[j+1]

            densInd = (df['density']>10**loN) & (df['density']<10**hiN)
            thisfit = np.zeros(len(header))
            thisfit[0] = a
            thisfit[1] = 1./a-1
            thisfit[2] = loN
            thisfit[3] = hiN
        
            cloud = df[densInd & spacInd & tempInd]
            thisfit[4] = len(cloud) / fullCloudSize

            cloud['speed'] = np.sqrt(cloud['vx']**2 + cloud['vy']**2 + cloud['vz']**2)

            if len(cloud)>1e4:
    
                # Radial distance and velocity
                cloud['vr'] = (cloud['x']*cloud['vx'] + cloud['y']*cloud['vy'] +
                                cloud['z']*cloud['vz'] ) / cloud['r']
                cloud['rMod'] = cloud['r']/rvir
                cloud['vrFrac'] = cloud['vr']/cloud['speed']
                cloud['vrFracAbs'] = np.abs(cloud['vr'])/cloud['speed']
            
                # Rho, the distance from rotated z-axis
                cloud['vrhoRot'] = np.sqrt(cloud['vxRot']**2 + cloud['vyRot']**2)

                # Same, but for rotated coordinates for sanity check
                cloud['rRot'] = np.sqrt(cloud['xRot']**2 + cloud['yRot']**2 + cloud['zRot']**2)
                cloud['vrRot'] = (cloud['xRot']*cloud['vxRot'] + cloud['yRot']*cloud['vyRot'] +
                                cloud['zRot']*cloud['vzRot'] ) / cloud['rRot']
                #cloud['vStat'] = cloud['vrhoRot'] / cloud['vzRot'] 

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
                cloud['valong'] = np.sqrt(cloud['along0']**2 +
                                cloud['along1']**2 + cloud['along2']**2)
                cloud['perp0'] = cloud['vx'] - cloud['along0']
                cloud['perp1'] = cloud['vy'] - cloud['along1']
                cloud['perp2'] = cloud['vz'] - cloud['along2']
                cloud['vperp'] = np.sqrt(cloud['perp0']**2 +
                                cloud['perp1']**2 + cloud['perp2']**2)
                
                # Get spherical coordinates
                cloud['phiRot'] = np.degrees(np.arctan(cloud['yRot']/cloud['xRot']))
                cloud['thetaRot'] = np.degrees(np.arccos(cloud['zRot']/cloud['rRot']))

                cloud['thetadot'] = ( (-1/np.sqrt( 1 - (cloud['zRot']/cloud['rRot'])**2)) * 
                                      (cloud['vzRot']/cloud['rRot'] - cloud['zRot']*cloud['vrRot']/cloud['rRot']**2))
                cloud['vthetaRot'] = cloud['rRot']*cloud['thetadot']
                cloud['phidot'] = ( (cloud['xRot']*cloud['vyRot'] - cloud['yRot']*cloud['vxRot']) / 
                                    (cloud['xRot']**2 + cloud['yRot']**2) )
                cloud['vphiRot'] = cloud['rRot']*np.sin(np.radians(cloud['thetaRot']))*cloud['phidot']


                # Mass and thermal properties 
                # Mass of the cell in solar Masses
                cloud['mass'] = (cloud['density']*u2g*pc2cm**3/g2M)*cloud['cell_size']**3  
                cloud['pressure'] = boltz*cloud['density']*cloud['temperature']
                cloud['thermalV'] = np.sqrt(8*boltz*cloud['temperature']/(u2g*np.pi)) * cm2km

                # Fill output array
                thisfit[5] = locM['dist'].std()
                thisfit[6] = param[0]
                thisfit[7] = param[1]
                thisfit[8] = fwhm
                thisfit[9] = st.rayleigh.std(loc=param[0],scale=param[1])

                index = 10
                for field in fields:
                    thisfit[index] = cloud[field].mean()
                    thisfit[index+1] = cloud[field].std()
                    thisfit[index+2] = thisfit[index+1]/thisfit[index]
                    thisfit[index+3] = np.average(cloud[field],weights=cloud['mass'])
                    thisfit[index+4] = cloud[field].median()
                    index += 5

            else:
                thisfit[5:] = np.NAN
                thisfit[10] = cloud['speed'].std()
                
            fit = np.vstack((fit,thisfit))
    

    fit = np.delete(fit,(0),axis=0)
    df = pd.DataFrame(fit,columns=header)
    if excludeCGM:
        outfile = 'simple_density_cuts_parameters_{0:s}_noCGM.h5'.format(tempLabels[i])
    else:
        outfile = 'simple_density_cuts_parameters_{0:s}.h5'.format(tempLabels[i])
    df.to_hdf(outfile,'data',mode='w')

                

















