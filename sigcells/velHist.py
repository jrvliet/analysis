'''
Plots phase space of LOS position and LOS velocity
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.stats.mstats import gmean
import sys

def plotHist(hs, xeds, yeds, ions, filename, xname, yname, cname ):

    mask = 0

    # Plot
    fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(1,4,figsize=(20.8, 5.2))
    axes = [ax1, ax2, ax3, ax4]


    for i in range(0,len(axes)):
        ion = ions[i]
        ax = axes[i]
        H = hs[i]
        if mask==1:
            H = np.ma.masked_where(H==0,H)
        else:
            minval = (np.ma.masked_where(H==0,H)).min()
            H[H<minval] = minval
    
        xedges = xeds[i]
        yedges = yeds[i]

        cm = mp.cm.get_cmap('viridis')
        mesh = ax.pcolormesh(xedges,yedges,H, cmap=cm)
        ax.set_xlabel(xname)
        ax.set_ylabel(yname)
        ax.set_title(ion)
        if 'normed' in filename:
            ax.set_ylim([0,1])

        if 'abs' in xname.lower():
            ax.set_xlim((0.0,xedges[-1]))
        else:
            ax.set_xlim((xedges[0],xedges[-1]))
    
        if 'abs' in yname.lower():
            ax.set_ylim((0.0,yedges[-1]))
        else: 
            ax.set_ylim((yedges[0],yedges[-1]))

        cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
        cbar.ax.get_yaxis().labelpad = 20
        cbar.ax.set_ylabel(cname, rotation=270, fontsize=12)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    fig.savefig(filename.replace('pdf','png'), bbox_inches='tight')






def binning_z(x, y, z, nbins, xlims, ylims, mode='hist'):
    '''
    Creates a 2D histogram of the data based on x and y
    Colored by z. How z is calcualted depends on mode:
        hist = z is ignored, bin values are number of cells
        mean = bin value is mean of z for cells in that bin
        gmean = bin value is geometric mean of z for cells in that bin
        std = standard deviation
        spread = difference between min and max
    '''

    # Make sure mode is valid
    validModes = ['hist', 'mean', 'gmean', 'std', 'spread']
    if mode not in validModes:
        print 'Nonvalid mode: {0:s}\nExiting...\n'.format(mode)
        sys.exit()

    # Make the bins
    xbins = np.linspace(xlims[0], xlims[1], nbins+1)
    ybins = np.linspace(ylims[0], ylims[1], nbins+1)

    # Determine what cells go in what bins
    xdig = np.digitize(x, xbins)
    ydig = np.digitize(y, ybins)

    # Fix the edge effects
    maxBinNum = len(xbins)
    for i in range(len(xdig)):
        if xdig[i]==maxBinNum:
            xdig[i] -= 1
        if ydig[i]==maxBinNum:
            ydig[i] -= 1

    # Create empty array
    h = np.zeros((nbins, nbins))

    # Loop through array
    for i in range(nbins):
        for j in range(nbins):
            # Find the indicies where x and y belong to this bin
            bits = np.bitwise_and( xdig==i+1, ydig==j+1 )
    
            if np.any(bits):
                if mode=='hist':
                    h[i,j] = np.sum( z[bits] )
                elif mode=='mean':
                    h[i,j] = np.mean( z[bits] )
                elif mode=='gmean':
                    h[i,j] = gmean( z[bits] )
                elif mode=='std':
                    h[i,j] = np.std( z[bits] )
                elif mode=='spread':
                    h[i,j] =  z[bits].max() - z[bits].min()
                else:
                    print 'Nonvalid mode: {0:s}\nExiting...\n'.format(mode)
                    sys.exit()

    # Flip, rotate, and log the array
    h = np.rot90(h)
    h = np.flipud(h)
    h = np.log10(h)
    h[h==-np.inf] = 0
    return h, xbins, ybins




ions = ['HI', 'MgII', 'CIV', 'OVI']

# Set binning parameters
numbins = 50
smin, smax = -220, 220
vmin, vmax = -250, 250
rmin, rmax = 0.0, 300
spmin, spmax = 0.0, 300
slims = (smin, smax)
vlims = (vmin, vmax)
absvlims = (0.0, vmax)
rlims = (rmin, rmax)
splims = (spmin, spmax)


printrange = 1

# Initialize arrays to hold binning
h_losv_s, xed_losv_s, yed_losv_s = [], [], []
h_abslosv_s, xed_abslosv_s, yed_abslosv_s = [], [], []
h_perpv_s, xed_perpv_s, yed_perpv_s = [], [], []
h_r_speed, xed_r_speed, yed_r_speed = [], [], []
h_r_vr, xed_r_vr, yed_r_vr = [], [], []

histos, xed, yed = [], [], []
histosabs, xedabs, yedabs = [], [], []
hs, xeds, yeds = [], [], []
h_rs, xed_rs, yed_rs = [], [], []
hperp, xedperp, yedperp = [], [], []
ss, vs = [], []
for ion in ions:
    
    print ion
    # Read in data    
    filename = '{0:s}_vlos.dat'.format(ion)
    l, s, v, n, t, Z, size, nIon, speed, vperp, r, vrad, col = np.loadtxt(filename, 
                        skiprows=1, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12), unpack=True)
    
    if printrange==1:
        print '\nFor {0:s}'.format(ion)
        print 'Range of LOS Length: {0:f}-{1:f}'.format(l.min(), l.max())
        print 'Range of S:          {0:f}-{1:f}'.format(s.min(), s.max())
        print 'Range of Vlos:       {0:f}-{1:f}'.format(v.min(), v.max())
        print 'Range of Density:    {0:f}-{1:f}'.format(n.min(), n.max())
        print 'Range of Temp:       {0:f}-{1:f}'.format(t.min(), t.max())
        print 'Range of AlphaZ:     {0:f}-{1:f}'.format(Z.min(), Z.max())
        print 'Range of Cellsize:   {0:f}-{1:f}'.format(size.min(), size.max())
        print 'Range of nIon:       {0:f}-{1:f}'.format(nIon.min(), nIon.max())
        print 'Range of Speed:      {0:f}-{1:f}'.format(speed.min(), speed.max())
        print 'Range of Vperp:      {0:f}-{1:f}'.format(vperp.min(), vperp.max())
        print 'Range of Distance:   {0:f}-{1:f}'.format(r.min(), r.max())
        print 'Range of Vrad:       {0:f}-{1:f}'.format(vrad.min(), vrad.max())
        print 'Range of Col Dense:  {0:f}-{1:f}'.format(col.min(), col.max())

    diff = [s[i] - v[i] for i in range(len(s))]
    ss.append(s)
    vs.append(v)

    # Scale LOS position so 0 is at the middle   
    mid = max(l) / 2.0
    for i in range(0,len(s)):
        s[i] = s[i]*l[i] - mid
        
    # Unlog all logged data points
    n    = np.array([ 10**i for i in n])
    t    = np.array([ 10**i for i in t])
    Z    = np.array([ 10**i for i in Z])
    nIon = np.array([ 10**i for i in nIon])
    col  = np.array([ 10**i for i in col])

    # Make histogram
    # Bin up the LOS vel and LOS location, colored by geo mean of column density
    h_losv_s_type = 'gmean'
    H, xedges, yedges = binning_z(v, s, col, numbins, vlims, slims, 'hist')
    xed_losv_s.append(xedges)
    yed_losv_s.append(yedges)
    h_losv_s.append(H)

    # Same as the one above, but use absolute value of LOS vel
    v = np.absolute(v)
    h_abslosv_s
    H, xedges, yedges = binning_z(v, s, col, numbins, absvlims, slims, 'hist')
    xed_abslosv_s.append(xedges)
    yed_abslosv_s.append(yedges)
    h_abslosv_s.append(H)

    # Bin the perpidcular component of velocity with LOS position
    h, xedges, yedges = binning_z(vperp, s, col, numbins, absvlims, slims, 'hist')
    xed_perpv_s.append(xedges)
    yed_perpv_s.append(yedges)
    h_perpv_s.append(h)

    # Bin the galactocentric distance and the speed
    h, xedges, yedges = binning_z(r, speed, col, numbins, rlims, splims, 'hist')
    xed_r_speed.append(xedges)
    yed_r_speed.append(yedges)
    h_r_speed.append(h)

    # Bin the radial velocity with the galactocentric distance
    h, xedges, yedges = binning_z(r, vrad, col, numbins, rlims, splims, 'hist')
    xed_r_vr.append(xedges)
    yed_r_vr.append(yedges)
    h_r_vr.append(h)
        


# Labels
losvel = 'LOS Velocity [km/s]'
abslosvel = 'Abs LOS Velocity [km/s]'
lospos = 'LOS Position [kpc]'
gColDense = 'Geo Mean Col Dense'
speed = 'Speed [km/s]'
hist = 'Number of cells'
dist = 'Galactocentric Distance [kpc]'
perpVel = 'Perpidicular Velocity [km/s]'
vradial = 'Radial Velocity [km/s]'

plotHist(h_losv_s, xed_losv_s, yed_losv_s, ions, 'vel2dHist_losv_s_hist.pdf', losvel, lospos, hist)
plotHist(h_abslosv_s, xed_abslosv_s, yed_abslosv_s, ions, 'vel2dHist_abslosv_s_hist.pdf', abslosvel, lospos, hist)
#plotHist(h, xeds, yeds, ions, 'vel2dHist_gmeanSpeed_new.pdf', losvel, speed, gColDense)
plotHist(h_perpv_s, xed_perpv_s, yed_perpv_s, ions, 'vel2dHist_vperp_s_hist.pdf', perpVel, losvel, hist)
plotHist(h_r_speed, xed_r_speed, yed_r_speed, ions, 'vel2dHist_r_s_speed_hist.pdf', dist, speed,  hist)
plotHist(h_r_vr, xed_r_vr, yed_r_vr, ions, 'vel2dHist_r_vr_hist.pdf', dist, vradial, hist)


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################



def make_histos(v, s, nIon, size, nbins, ion):

    '''
    Bins up the x and y data into nbins
    Value is each bin is the geometric mean of column density
    '''
    
    # Find column density
    column = np.zeros(len(size))
    for i, (n,l) in enumerate(zip(nIon, size)):
        # Take the cube root of the cell length and convert from kpc to cm
        length = l**(1./3.) * 3.086e21
        column[i] = (10**n) * length
    
    print min(column)
    print max(column)
    vmin = -100
    vmax = 100
    smin = -220
    smax = 220  

    print 'Making histogram'
    H, xed, yed = np.histogram2d(v,s,bins=nbins, range=[[vmin,vmax],[smin,smax]])
    h = np.zeros_like(H)
    print 'Histogram done\n'

    for i in range(0,H.shape[0]):
        for j in range(0,H.shape[1]):
            # #rows = shape[0]
            # #cols = shape[1]
            vmin = xed[j]
            vmax = xed[j+1]
            smin = yed[i]
            smax = yed[i+1]
            val = []
            for k, (vel,sloc) in enumerate(zip(v,s)):
                if vel>vmin and vel<vmax and sloc>smin and sloc<smax:
                    val.append(column[k])
            print np.log10(min(val))
            print np.log10(max(val))
            print np.log10(np.mean(val))
            print np.log10(gmean(val))
            print vmin, vmax
            print smin, smax
            h[i,j] = gmean(val)
            print h[i,j]
    h = np.log10(h)
    np.savetxt('{0:s}_velHist.out'.format(ion), h)
    return h,xed,yed


def hist(x, y, numbins, ion):

    # Make the bins
    xbins = np.linspace(xlims[0], xlims[1], nbins+1)
    ybins = np.linspace(ylims[0], ylims[1], nbins+1)

    # Determine what cells go in what bins
    xdig = np.digitize(x, xbins)
    ydig = np.digitize(y, ybins)

    # Fix the edge effects
    maxBinNum = len(xbins)
    for i in range(len(xdig)):
        if xdig[i]==maxBinNum:
            xdig[i] -= 1
        if ydig[i]==maxBinNum:
            ydig[i] -= 1
    
    # Create empty array
    h = np.zeros((nbins, nbins))

    # Loop through array
    for i in range(nbins):
        for j in range(nbins):
            # Find the indicies where x and y belong to this bin
            bits = np.bitwise_and( xdig==i+1, ydig==j+1)
            if True in bits:
                h[i,j] = np.log10( np.sum( bits ) )

    h = np.rot90(h)
    h = np.flipud(h)
    return h, xbins, ybins












def gmean_speed_bin(x, y, nbins,xlims,ylims,ion,speed):
    '''
    Bins up the data according to x and y
    Value in each bin is the geometric mean of the speed
    '''

    x = speed
    # Make the bins
    xbins = np.linspace(xlims[0], xlims[1], nbins+1)
    ybins = np.linspace(ylims[0], ylims[1], nbins+1)

    # Determine what cells go in what bins
    xdig = np.digitize(x, xbins)
    ydig = np.digitize(y, ybins)

    # Fix the edge effects
    maxBinNum = len(xbins)
    for i in range(len(xdig)):
        if xdig[i]==maxBinNum:
            xdig[i] -= 1
        if ydig[i]==maxBinNum:
            ydig[i] -= 1
    
    # Create empty array
    h = np.zeros((nbins, nbins))

    # Loop through array
    for i in range(nbins):
        for j in range(nbins):
            # Find the indicies where x and y belong to this bin
            bits = np.bitwise_and( xdig==i+1, ydig==j+1)
            if True in bits:
                h[i,j] = np.log10( np.sum( speed[bits] ) )

    h = np.rot90(h)
    h = np.flipud(h)
    return h, xbins, ybins




def gmean_bin(x, y, nIon, size, nbins, xlims, ylims, ion):
    '''
    Bins up the data according to x and y
    Value in each bin is the geometric mean of the column
    density contribution of cells in that bin, as determined
    by multiplying nIon by size**1/3
    '''

    # Calculate the column density
    column = np.zeros(len(size))
    f = open('{0:s}_column.out'.format(ion), 'w')
    for i, (n,l) in enumerate(zip(nIon, size)):
        # Take the cube root of the cell length and convert from kpc to cm
        length = l**(1./3.) * 3.086e21
        col = n * length
        column[i] = col
        f.write('{0:.4e}\n'.format(col))
    f.close()        

    # Make the bins
    xbins = np.linspace(xlims[0], xlims[1], nbins+1)
    ybins = np.linspace(ylims[0], ylims[1], nbins+1)

    # Determine what cells go in what bins
    xdig = np.digitize(x, xbins)
    ydig = np.digitize(y, ybins)

    # Fix the edge effects
    maxBinNum = len(xbins)
    for i in range(len(xdig)):
        if xdig[i]==maxBinNum:
            xdig[i] -= 1
        if ydig[i]==maxBinNum:
            ydig[i] -= 1
    
    # Create empty array
    h = np.zeros((nbins, nbins))

    # Loop through array
    for i in range(nbins):
        for j in range(nbins):
            # Find the indicies where x and y belong to this bin
            bits = np.bitwise_and( xdig==i+1, ydig==j+1)
            if True in bits:
                h[i,j] = np.log10( gmean( column[bits] ) )

    h = np.rot90(h)
    h = np.flipud(h)
    np.savetxt('{0:s}_velHist.out'.format(ion), h)
    print 'Max of h: ', np.max(h)
    print 'Mean of h: ', np.mean(h)
    return h, xbins, ybins
