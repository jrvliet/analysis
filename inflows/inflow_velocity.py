

'''
Determines the coherence of the inflows
Focus on the poster child inflow in vela2b-27
at a=0.490
'''


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.linalg as sl
import scipy.stats as st
import sys

def mkHist(ax,x,lims,numbins,label,color,xlab, log=False):

    log = True
    l = ax.hist(x,bins=numbins,range=lims,histtype='step',
                log=log,label=label,color=color)
    ax.set_xlabel(xlab)
    ax.set_ylabel('Counts')

    return l


pd.options.mode.chained_assignment = None

plotting = 0

loT, hiT = 10**3.25, 10**4.5

minN, maxN, numbins = -6, -2, 7
nBins = np.linspace(minN,maxN,numbins)

maxR = 6 
rBins = np.arange(0,maxR+1,1)

expns = [0.490]
expns = np.arange(0.200,0.500,0.010)

baseloc = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/'
baseloc = '/home/jacob/research/velas/vela2b/vela27/a{0:.3f}/'
basename = '{0:s}vela2b-27_GZa{1:.3f}.h5'
baserot = '{0:s}rotmat_a{1:.3f}.txt'


vxLims, vyLims, vzLims = [],[],[]
vrLims, vthetaLims, vphiLims = [],[],[]
vspeedLims, valongLims, vperpLims = [],[],[]
distLims = []


header = ['a','z','sigma_valong','sigma_vperp','sigma_r','sigma_theta','sigma_phi']
vel_nbin1_rbinall = np.zeros(len(header))
vel_nbin2_rbinall = np.zeros(len(header))
vel_nbin3_rbinall = np.zeros(len(header))
vel_nbin4_rbinall = np.zeros(len(header))
vel_nbin5_rbinall = np.zeros(len(header))
vel_nbin6_rbinall = np.zeros(len(header))
vel_nbins = [vel_nbin6_rbinall,vel_nbin6_rbinall,vel_nbin6_rbinall,vel_nbin6_rbinall,vel_nbin6_rbinall,vel_nbin6_rbinall]

for a in expns:

    print('a = {0:.3f}'.format(a))
    redshift = 1./a - 1
    
    loc = baseloc.format(a)
    fname = basename.format(loc,a)
    frot = baserot.format(loc,a)

    with open(frot, 'r') as f:
        f.readline()
        rvir = float(f.readline().split()[3])
    
    df = pd.read_hdf(fname, 'data')

    # Loop over density bins
    for i in range(numbins-1):

        print('\tDensity range = {0:.3f} - {1:.3f}'.format(nBins[i],nBins[i+1]))
        
        loN = 10**nBins[i]
        hiN = 10**nBins[i+1]

        index = ( (df['temperature']>loT) & (df['temperature']<hiT) & 
                    (df['density']>loN) & (df['density']<hiN) &
                    (df['x']>0) & (df['z']>0) & (np.abs(df['y'])<300) )

        cloud = df[index]
        cloudLoc = cloud[['x','y','z']]
        print(len(cloudLoc))

        if len(cloudLoc)>0:
            #####################################################
            # Start working with velocities
            cloud['speed'] = np.sqrt(cloud['vx']**2 + cloud['vy']**2 + cloud['vz']**2 )

            # Determine the speed of each cell
            loSpeed = cloud['speed'].min()
            hiSpeed = cloud['speed'].max()

            # Calculate the spherical coordinates
            cloud['r'] = np.sqrt(cloud['x']**2 + cloud['y']**2 + cloud['z']**2 )
            cloud['theta'] = np.arctan2(cloud['y'],cloud['x'])
            cloud['phi'] = np.arccos(cloud['z']/cloud['r'])
            
            # Calculate the spherical velocities
            cloud['vr'] = (cloud['x']*cloud['vx'] + cloud['y']*cloud['vy'] + 
                            cloud['z']*cloud['vz']) / cloud['r']

            cloud['thetadot'] = ( (cloud['x']*cloud['vy'] - cloud['y']*cloud['vx']) /
                                (cloud['x']**2 + cloud['y']**2 ) )
            cloud['vtheta'] = cloud['r']*np.sin(cloud['phi'])*cloud['thetadot']

            cloud['phidot'] = ( (-1 / np.sqrt( 1 - (cloud['z']/cloud['r'])**2)) * 
                                (cloud['vz']/cloud['r'] - cloud['z']*cloud['vr']/cloud['r']**2) )
            cloud['vphi'] = cloud['r']*cloud['phidot']

            distLims.append(cloud['r'].max()/rvir)


            #######################################################

            # Perform a PCA to find the line of best fit
            # Start by normalizing the cells
            locM = cloudLoc - cloudLoc.mean()

            # Get the covariance matrix
            covar = locM.cov()

            # Determine the single value decomposition of the covariance matrix
            (u,s,v) = sl.svd(covar)

            # The first column of u is the directional vector of the line
            # of best fit for the cloud
            u0 = u[0,0]
            u1 = u[1,0]
            u2 = u[2,0]

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


            # Fill the arrays
            vel = np.zeros(len(header))
            vel[0] = a
            vel[1] = 1./a - 1
            vel[2] = cloud['along'].std()
            vel[3] = cloud['perp'].std()
            vel[4] = cloud['vr'].std()
            vel[5] = cloud['vtheta'].std()
            vel[6] = cloud['vphi'].std()
            vel_nbins[i] = np.vstack((vel_nbins[i],vel))

            #######################################################
        
            if plotting==1:
                # Plot historgram of dist
                fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(15,15))
                axes = (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9)
                log = False

                colors = ['blue','red','green','black','purple','orange']
                labels = ['0<r<1','1<r<2','2<r<3','3<r<4','4<r<5','5<r<6']

                for i in range(len(rbins)-1):

                    rmin = rbins[i]*rvir
                    rmax = rbins[i+1]*rvir

                    clInds = ((cloud['r']>=rmin) & (cloud['r']<rmax))
                    cl = cloud[clInds]

                    # Plot cartesian velocities
                    lims = [cl['vx'].min(), cl['vx'].max()]
                    vxLims.append(lims[0])
                    vxLims.append(lims[1])
                    l1 = mkHist(ax1, cl['vx'], lims, numbins, labels[i], colors[i], log, )

                    lims = [cl['vy'].min(), cl['vy'].max()]
                    vyLims.append(lims[0])
                    vyLims.append(lims[1])
                    l2 = mkHist(ax2, cl['vy'], lims, numbins, labels[i], colors[i], 'vy', log)

                    lims = [cl['vz'].min(), cl['vz'].max()]
                    vzLims.append(lims[0])
                    vzLims.append(lims[1])
                    l3 = mkHist(ax3, cl['vz'], lims, numbins, labels[i], colors[i],'vz', log)


                    # Plot spherical velocities
                    lims = [cl['vr'].min(), cl['vr'].max()]
                    vrLims.append(lims[0])
                    vrLims.append(lims[1])
                    l4 = mkHist(ax4, cl['vr'], lims, numbins, labels[i], colors[i], 'vr', log)
                    #mkHist(ax4, df['vr'], lims, numbins, 'all', 'vr', log)

                    lims = [cl['vphi'].min(), cl['vphi'].max()]
                    vphiLims.append(lims[0])
                    vphiLims.append(lims[1])
                    l5 = mkHist(ax5, cl['vphi'], lims, numbins, labels[i], colors[i],'vphi', log)
                    #mkHist(ax5, df['vphi'], lims, numbins, 'all', 'vphi', log)

                    lims = [cl['vtheta'].min(), cl['vtheta'].max()]
                    vthetaLims.append(lims[0])
                    vthetaLims.append(lims[1])
                    l6 = mkHist(ax6, cl['vtheta'], lims, numbins, labels[i], colors[i], 'vtheta', log)
                    #mkHist(ax6, df['vtheta'], lims, numbins, 'all', 'vtheta', log)
                    
                    
                    # Plot the speeds
                    log = False
                    lims = [cl['speed'].min(), cl['speed'].max()]
                    vspeedLims.append(lims[0])
                    vspeedLims.append(lims[1])
                    l7 = ax7.scatter(cl['along'],cl['vr'])
                    ax7.set_xlabel('Valong')
                    ax7.set_ylabel('Vr')
                    ax7Lims = ax7.get_ylim()
                    #l7 = mkHist(ax7, cl['speed'], lims, numbins, labels[i], colors[i], 'Speed')

                    lims = [cl['along'].min(), cl['along'].max()]
                    valongLims.append(lims[0])
                    valongLims.append(lims[1])
                    l8 = mkHist(ax8, cl['along'], lims, numbins, labels[i], colors[i], 'Along')

                    lims = [cl['perp'].min(), cl['perp'].max()]
                    vperpLims.append(lims[0])
                    vperpLims.append(lims[1])
                    l9 = mkHist(ax9, cl['perp'], lims, numbins, labels[i], colors[i], 'V Perp')
                    
                
                # Set limits
                for ax in axes:
                    ax.set_ylim(ymin=0)
                ax7.set_ylim(ax7Lims)

                ax1.set_xlim([-400,400])
                ax2.set_xlim([-400,400])
                ax3.set_xlim([-400,400])

                ax4.set_xlim([-400,400])
                ax5.set_xlim([-400,400])
                ax6.set_xlim([-400,400])
                
                ax7.set_xlim([0,500])
                ax8.set_xlim([0,500])
                ax9.set_xlim([0,500])

                

                lines = (l1,l2,l3,l4,l5,l6,l7,l8,l9)
                #fig.legend(lines,labels,'upper center')

                # Save
                fig.suptitle('z = {0:.3f}'.format(redshift), fontsize=26)
                fig.tight_layout()
                fig.subplots_adjust(top=0.95)
                s = 'vela2b-27_a{0:.3f}_inflowVelocity.png'.format(a)
                fig.savefig(s, bbox_inches='tight', dpi=300)


        else:
            # Fill the arrays
            vel = np.zeros(len(header))
            vel_nbins[i] = np.vstack((vel_nbins[i],vel))
    
    #######################################################

outfilename = 'velocity_dispersions.h5'
store = pd.HDFStore(outfilename)
for i in range(len(vel_nbins)):
    vel = np.delete(vel_nbins[i],(0),axis=0)
    df = pd.DataFrame(vel,columns=header)
    name = 'nbin{0:d}'.format(i+1)
    store[name] = df
store.close()
    
    
if plotting==1:
    print('Vx     = {0:.3f} - {1:.3f}'.format(min(vxLims), max(vxLims)))
    print('Vy     = {0:.3f} - {1:.3f}'.format(min(vyLims), max(vyLims)))
    print('Vz     = {0:.3f} - {1:.3f}'.format(min(vzLims), max(vzLims)))
    print('Vr     = {0:.3f} - {1:.3f}'.format(min(vrLims), max(vrLims)))
    print('Vphi   = {0:.3f} - {1:.3f}'.format(min(vphiLims), max(vphiLims)))
    print('Vtheta = {0:.3f} - {1:.3f}'.format(min(vthetaLims), max(vthetaLims)))
    print('Vspeed = {0:.3f} - {1:.3f}'.format(min(vspeedLims), max(vspeedLims)))
    print('Valong = {0:.3f} - {1:.3f}'.format(min(valongLims), max(valongLims)))
    print('Vperp  = {0:.3f} - {1:.3f}'.format(min(vperpLims), max(vperpLims)))


    print('Max Distance = {0:.3f}'.format(max(distLims)))














