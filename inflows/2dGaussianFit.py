
from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
%matplotlib inline
import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.ticker as ticker
pd.options.mode.chained_assignment = None
import scipy.optimize as opt



def fmt(x,pos):
    return '{0:.2f}'.format(x)

def mkHist(x,y,z,stat,xlabel,ylabel):
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
 
    return h,xedges,yedges,binnumber


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2.)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y



def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = opt.leastsq(errorfunction, params)
    return p







kpc2km = 3.086e16 
pc2cm = 3.086e18 
mH = 1.6737e-24 
mSun = 1.989e33 
s2yr = 3.154e7

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'
expn = np.arange(0.200,0.500,0.01)

for a in expn:
print(a)
rname = dataloc+rotname.format(a)
with open(rname) as f:
    f.readline()
    rvir = float(f.readline().split()[3])

fname = dataloc + filename.format(a)
df = pd.read_hdf(fname,'data')
df['mass'] = df['density']*mH*(df['cell_size']*pc2cm)**3 / mSun

loT,hiT = 3.5,4.5
loN,hiN = -5.5,-2.5
tempInds = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
denseInds = (df['density']>10**loN) & (df['density']<10**hiN)
spaceInds = (df['theta']<80) & (df['r']>0.5*rvir)
fil = df[tempInds & denseInds & spaceInds]

results = mkHist(fil['xRot']/rvir,fil['yRot']/rvir,fil['mass'],'count','xRot','yRot')



h = results[0]
x = results[1][:-1]
y = results[2][:-1]



h[h.mask] = 0


params = fitgaussian(h)


binrange = [[-3,3],[-3,3]]

x0 = params[1]
y0 = params[2]
xsigma = params[3]
ysigma = params[4]












