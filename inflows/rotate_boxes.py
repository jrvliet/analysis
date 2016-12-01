
'''
Creates rotated gasboxes
'''

from __future__ import print_function
import numpy as np
import pandas as pd


dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'

rotmat = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

rot = np.loadtxt('inflow_rot_mat.dat')

expns = np.arange(0.200,0.500,0.01)
expns = np.arange(0.490,0.550,0.01)

for a in expns:
    
    print(a)

    fname = dataloc+filename.format(a)
    rname = dataloc+rotmat.format(a)

    with open(rname) as f:
        f.readline()
        l = f.readline().split()
        a = np.array(map(float,l[4:-3])).reshape(3,3)

    cl = pd.read_hdf(fname, 'data')
    cl['r'] = np.sqrt(cl['x']**2 + cl['y']**2 + cl['z']**2)
    cl['theta'] = np.degrees(np.arccos(cl['z']/cl['r']))
    cl['phi'] = np.degrees(np.arctan2(cl['y'],cl['x']))
    cl['vr'] = (cl['x']*cl['vx'] + cl['y']*cl['vy'] + cl['z']*cl['vz']) / cl['r'] 

    # Start with the inflow rotation
    # Rotate positions
    clRot = cl[['x','y','z']].dot(rot) 
    clRot.rename( columns={0:'xRot', 1:'yRot', 2:'zRot'}, inplace=True)

    clRot['rRot'] = np.sqrt(clRot['xRot']**2 + clRot['yRot']**2 + clRot['zRot']**2)
    clRot['thetaRot'] = np.degrees(np.arccos(clRot['zRot']/clRot['rRot']))
    clRot['phiRot'] = np.degrees(np.arctan2(clRot['yRot'],clRot['xRot']))

    cl = cl.join(clRot)

    # Rotate velocities
    clRot = cl[['vx','vy','vz']].dot(rot)
    clRot.rename( columns={0:'vxRot', 1:'vyRot', 2:'vzRot'}, inplace=True)
    cl = cl.join(clRot)

#
#    cl['vrRot'] = (cl['xRot']*cl['vxRot'] + cl['yRot']*cl['vyRot'] + 
#                    cl['zRot']*cl['vzRot']) / cl['rRot'] 
#
#    # Now do the galaxy rotation
#    clRot = cl[['x','y','z']].dot(a)
#    clRot.rename( columns={0:'xGal', 1:'yGal', 2:'zGal'}, inplace=True)
#    
#    clRot['rGal'] = np.sqrt(clRot['xGal']**2 + clRot['yGal']**2 + clRot['zGal']**2)
#    clRot['thetaGal'] = np.degrees(np.arccos(clRot['zGal']/clRot['rGal']))
#    clRot['phiGal'] = np.degrees(np.arctan2(clRot['yGal'],clRot['xGal']))
#
#    cl = cl.join(clRot)
#
#    # Rotate velocities
#    clRot = cl[['vx','vy','vz']].dot(a)
#    clRot.rename( columns={0:'vxGal', 1:'vyGal', 2:'vzGal'}, inplace=True)
#
#    cl = cl.join(clRot)
#    cl['vrGal'] = (cl['xGal']*cl['vxGal'] + cl['yGal']*cl['vyGal'] + 
#                    cl['zGal']*cl['vzGal']) / cl['rGal'] 
#
#
#    # Get the spherical velocities in each rotation
#    suffix = ['','Rot','Gal']
#    for suf in suffix:
#        
#        x = cl['x'+suf]
#        y = cl['y'+suf]
#        z = cl['z'+suf]
#        r = cl['r'+suf]
#        vx = cl['vx'+suf]
#        vy = cl['vy'+suf]
#        vz = cl['vz'+suf]
#        vr = cl['vr'+suf]
#        
#        theta = cl['theta'+suf]
#
#        cl['thetadot'] = ( (-1 / np.sqrt( 1 - (z/r)**2)) *
#                           (vz/r - z*vr/r**2) )
#        cl['phidot'] = (x*vy - y*vx) / (x**2 + y**2) 
#        cl['vtheta'+suf] = r*cl['thetadot']
#        cl['vphi'+suf] = r*np.sin(np.radians(theta))*cl['phidot']
#     

    #outname = fname.replace('.h5','.fullrot.h5')
    outname = fname.replace('.h5','.rot.h5')
    cl.to_hdf(outname,'data',mode='w')

