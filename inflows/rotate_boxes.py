
'''
Creates rotated gasboxes
'''

from __future__ import print_function
import numpy as np
import pandas as pd


dataloc = '/mnt/cluster/abs/cgm/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.h5'

rot = np.loadtxt('inflow_rot_mat.dat')

expns = np.arange(0.200,0.500,0.01)

for a in expns:
    
    print(a)

    fname = dataloc+filename.format(a)

    cl = pd.read_hdf(fname, 'data')

    clRot = cl[['x','y','z']].dot(rot) 
    clRot.rename( columns={0:'xRot', 1:'yRot', 2:'zRot'}, inplace=True)

    clRot['r'] = np.sqrt(clRot['xRot']**2 + clRot['yRot']**2 + clRot['zRot']**2)
    clRot['theta'] = np.degrees(np.arccos(clRot['zRot']/clRot['r']))
    clRot['phi'] = np.degrees(np.arctan2(clRot['yRot'],clRot['xRot']))

    cl = cl.join(clRot)

    outname = fname.replace('.h5','.rot.h5')
    cl.to_hdf(outname,'data',mode='w')

