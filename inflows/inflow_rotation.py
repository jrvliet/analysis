
'''
Rotates inflow to z-axis
'''


from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import scipy.linalg as sl
import scipy.stats as st
import matplotlib.pyplot as plt

def fit_line2d(x,y):
    slope, intercept, r, p, stderr = st.linregress(x,y)
    return slope, intercept

def fit_line3d(d):
    '''
    Fits a line to the 3d distribution of cells using PCA
    '''
    dloc = d[['x','y','z']]
    locM = dloc - dloc.mean()
    covar = locM.cov()
    (u,s,v) = sl.svd(covar)
    u0 = u[0,0]
    u1 = u[1,0]
    u2 = u[2,0]
    return [u0,u1,u2]

def vector_angle(a,b):
    '''
    Returns the angle between two vectors, a and b
    '''
    top = np.dot(a,b)
    lena = np.sqrt(np.dot(a,a))
    lenb = np.sqrt(np.dot(b,b))
    theta = np.arccos(top/(lena*lenb))
    return theta

def plotting(d,savename):
    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize(15,5))
    
    # Plot cells
    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
    ax1.scatter(d['x']/rvir,d['y']/rvir,color='b',marker='o',alpha=0.01)
    ax2.scatter(d['x']/rvir,d['z']/rvir,color='r',marker='o',alpha=0.01)
    ax3.scatter(d['y']/rvir,d['z']/rvir,color='g',marker='o',alpha=0.01)
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax2.set_xlabel('x')
    ax2.set_ylabel('z')
    ax3.set_xlabel('y')
    ax3.set_ylabel('z')
    
    ax1.set_xlim([-3,3])
    ax1.set_ylim([-3,3])
    ax2.set_xlim([-3,3])
    ax2.set_ylim([-3,3])
    ax3.set_xlim([-3,3])
    ax3.set_ylim([-3,3])
    fig.tight_layout()

    fig.savefig(savename,bbox_inches='tight',dpi=300)

def plotting_with_line(d,vector,center,savename):
    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
    ax1.scatter(d['x']/rvir,d['y']/rvir,color='b',marker='o',alpha=0.01)
    ax2.scatter(d['x']/rvir,d['z']/rvir,color='r',marker='o',alpha=0.01)
    ax3.scatter(d['y']/rvir,d['z']/rvir,color='g',marker='o',alpha=0.01)

    # Plot the line
    t = np.linspace(-1000,1000,10000)
    x, y, z = [],[],[]
    for i in t:
        x.append( (vector[0]*i + center[0]) / rvir)
        y.append( (vector[1]*i + center[1]) / rvir)
        z.append( (vector[2]*i + center[2]) / rvir)
        
    ax1.plot(x,y,'k-',linewidth=3)
    ax2.plot(x,z,'k-',linewidth=3)
    ax3.plot(y,z,'k-',linewidth=3)
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax2.set_xlabel('x')
    ax2.set_ylabel('z')
    ax3.set_xlabel('y')
    ax3.set_ylabel('z')
    
    ax1.set_xlim([-6,6])
    ax1.set_ylim([-6,6])
    ax2.set_xlim([-6,6])
    ax2.set_ylim([-6,6])
    ax3.set_xlim([-6,6])
    ax3.set_ylim([-6,6])
    fig.tight_layout()

    fig.savefig(savename,bbox_inches='tight',dpi=300)


def plot_gradients(fullCloud,fullBox,field,savename):

    # Plot the metallicity gradients 
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
    ax1.scatter(fullCloud['theta'],np.log10(fullCloud[field]),color='green',marker='o',alpha=0.01)
    ax2.scatter(fullBox['theta'],np.log10(fullBox[field]),color='plum',marker='o',alpha=0.01)

    # Fit lines
    t = np.linspace(0,180,1000)
    mCloud,bCloud = fit_line2d(fullCloud['theta'],np.log10(fullCloud[field]))
    mBox,bBox = fit_line2d(fullBox['theta'],np.log10(fullBox[field]))
    y = mCloud*t + bCloud
    ax1.plot(t,y,'k')
    y = mBox*t + bBox
    ax2.plot(t,y,'k')
    print('\nFits for {0:s}'.format(field))
    print('mCloud = {0:.3f}\tbCloud = {1:.3f}'.format(mCloud,bCloud))
    print('mBox   = {0:.3f}\tbBox   = {1:.3f}'.format(mCloud,bCloud))

    # Annotate fit parameters
    yloc = fullCloud[field].max()
    s = 'm={0:.3f}\nb={1:.3f}'.format(mCloud,bCloud)
    ax1.annotate(s=s,xy=[yloc,1])
    yloc = fullBox[field].max()
    s = 'm={0:.3f}\nb={1:.3f}'.format(mBox,bBox)
    ax2.annotate(s=s,xy=[yloc,1])
    

    # Polish plots
    ax1.set_xlabel('theta [deg]')
    ax2.set_xlabel('theta [deg]')
    ax1.set_ylabel(field)
    ax2.set_ylabel(field)
    ax1.set_title('Cloud')
    ax2.set_title('Full Box')

    fig.tight_layout()
    fig.savefig(savename, bbox_inches='tight', dpi=300)





dataloc = '/home/jacob/research/velas/vela2b/vela27/a{0:.3f}/'
dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/a{0:.3f}/'
filename = 'vela2b-27_GZa{0:.3f}.h5'
rotfile = dataloc+'rotmat_a{0:.3f}.txt'


# Define phase
loT, hiT = 10**3.5, 10**4.5
loN, hiN = 10**-4, 10**-3.75

# Read in rotation matrix
rot = np.loadtxt('inflow_rot_mat.dat')

expns = np.arange(0.200,0.500,0.01)
expns = np.arange(0.500,0.550,0.01)

for a in expns:

    print(a)
    with open(rotfile.format(a),'r') as f:
        f.readline()
        rvir = float(f.readline().split()[3])


    # Read in data
    fname = dataloc.format(a) + filename.format(a)
    df = pd.read_hdf(fname, 'data')

    # Select inflow
    tempInd = (df['temperature']<hiT) & (df['temperature']>loT)
    denseInd = (df['density']<hiN) & (df['density']>loN)
    spaceInd = (df['x']<0) & (df['z']>0) & (np.abs(df['y'])<300)

    cl = df[ tempInd & denseInd & spaceInd ]

    # Fit a line
    u = fit_line3d(cl)

    # Get the rotation angle
    # Will rotate around the y-axis to make the inflow lie along 
    # the positive z-axis. To do this the angle between the
    # line of best fit and the current x-axis is needed
    x = [1,0,0]
    theta = vector_angle(x,u)

    # This actually will put the inflow along the negative z-axis, 
    # so need to add pi to it
    theta += np.pi

    # Create the rotation matrix
    rot = np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-1*np.sin(theta),0,np.cos(theta)]])
    np.savetxt('inflow_rot_mat_a{0:.3f}.dat'.format(a),rot)

    # Rotate the coordinates
    clRot = cl[['x','y','z']].dot(rot)
    clRot.rename( columns={0:'xRot', 1:'yRot', 2:'zRot'}, inplace=True )

    # Get the spherical coordinates of the inflow
    clRot['r'] = np.sqrt(clRot['xRot']**2 + clRot['yRot']**2 + clRot['zRot']**2)
    clRot['theta'] = np.degrees(np.arccos(clRot['zRot']/clRot['r']))
    clRot['phi'] = np.degrees(np.arctan2(clRot['yRot'],clRot['xRot']))

    # Join this to the full dataset for the inflow 
    fullCloud = cl.join(clRot)

    # Rotate the entire box
    boxLocRot = df[['x','y','z']].dot(rot)
    boxLocRot.rename(columns={0:'xRot', 1:'yRot', 2:'zRot'}, inplace=True)

    # Get the spherical coordinates for the full, rotated box
    boxLocRot['r'] = np.sqrt(boxLocRot['xRot']**2 + boxLocRot['yRot']**2 + boxLocRot['zRot']**2)
    boxLocRot['theta'] = np.degrees(np.arccos(boxLocRot['zRot']/boxLocRot['r']))
    boxLocRot['phi'] = np.degrees(np.arctan2(boxLocRot['yRot'],boxLocRot['xRot']))

    # Join this to the full, rotated box
    fullBox = df.join(boxLocRot)

    # Write to file
    foutname = filename.replace('.h5','.rot.h5')
    foutname = fname.replace('.h5','.rot.h5')
    fullBox.to_hdf(foutname, 'data', mode='w')

    # Plot gradients
    #plot_gradients(fullCloud,fullBox,'SNII','inflow_SNII_thetaGradient.png')
    #plot_gradients(fullCloud,fullBox,'density','inflow_nH_thetaGradient.png')


















