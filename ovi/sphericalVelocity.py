
'''
Calculates the spherical velocities in the galaxy's coordinate 
system of all OVI absorbing cells
'''

from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def azimuthal(phi):
    if phi<90:
        return phi
    elif phi<180:
        return 180.-phi
    elif phi<270:
        return phi-180.
    else:
        return 360.-phi
    

def rotate(row,rot):
    '''
    Converts the coordinates to galaxy frame
    '''

    # Coordinates
    row['xGal'] = (row['x']*rot['a11'] + row['y']*rot['a12'] + 
                   row['z']*rot['a13']).values[0]
    row['yGal'] = (row['x']*rot['a21'] + row['y']*rot['a22'] + 
                   row['z']*rot['a23']).values[0]
    row['zGal'] = (row['x']*rot['a31'] + row['y']*rot['a32'] + 
                   row['z']*rot['a33']).values[0]
    
    # Velocities
    row['vxGal'] = (row['vx']*rot['a11'] + row['vy']*rot['a12'] + 
                   row['vz']*rot['a13']).values[0]
    row['vyGal'] = (row['vx']*rot['a21'] + row['vy']*rot['a22'] + 
                   row['vz']*rot['a23']).values[0]
    row['vzGal'] = (row['vx']*rot['a31'] + row['vy']*rot['a32'] + 
                   row['vz']*rot['a33']).values[0]
    return row


def spherical(row):
    '''
    Converts the cartesian velocities to spherical velocities
    '''

    # Coordinates
    row['r'] = np.sqrt(row['xGal']**2 + row['yGal']**2 + row['zGal']**2)
    row['theta'] = np.arctan2(row['yGal'],row['xGal'])
    row['phi'] = np.arccos(row['zGal']/row['r'])

    # Velocities
    row['vr'] = (row['xGal']*row['vxGal'] + row['yGal']*row['vyGal'] + 
                 row['zGal']*row['vzGal']) / row['r']
    row['thetadot'] = ((row['xGal']*row['vyGal'] - row['yGal']*row['vxGal']) / 
                      (row['xGal']**2 + row['yGal']**2))
    row['phidot'] = ( (-1/np.sqrt(1.-(row['zGal']/row['r'])**2)) * 
                    (row['vzGal']/row['r'] - row['zGal']*row['vr']/row['r']**2))
    row['vtheta'] = row['r']*np.sin(row['phi'])*row['thetadot']
    row['vphi'] = row['r']*row['phidot']
    
    row['phi'] = 180./np.pi*row['phi']
    row['theta'] = 180./np.pi*row['theta']

    return row

def cell_plane_cut(phi,planeAzCut):
    
    if phi>=90.-planeAzCut and phi<=90.+planeAzCut:
        return 1
    else:
        return 0

def cell_outflow_cut(phi,outflowAzCut):
    if phi<=outflowAzCut or phi>=180.-outflowAzCut:
        return 1
    else:
        return 0


def los_plane_cut(losnum,planeAzCut,lines):
    losPhi = lines['phi'].iloc[int(losnum)-1]
    losAz = azimuthal(losPhi)
    if losAz<=planeAzCut:
        return 1
    else:
        return 0

def los_outflow_cut(losnum,outflowAzCut,lines):
    losPhi = lines['phi'].iloc[int(losnum)-1]
    losAz = azimuthal(losPhi)
    if losAz>=outflowAzCut:
        return 1
    else:
        return 0

baseloc = '/mnt/cluster/abs/cgm/vela2b/'
subloc = 'vela{0:d}/a{1:.3f}/i90/OVI/'
cellFileName = 'vela2b-{0:d}.{1:.3f}.OVI.i90.abs_cells.h5'
linesName = 'lines.info'
rotFileName = '../../rotmat_a{0:.3f}.txt'
outFileName = 'vela2b-{0:d}.a{1:.3f}.i90_OVIcellsGalFrame.h5'


cols = 'LOS cellID x y z vx vy vz'.split()
linesHeader = 'losnum impact phi inclination'.split()
outHeader = ['cellID','x','y','z','vx','vy','vz',
             'r','theta','phi','vr','vtheta','vphi',
             'cellPlane','cellOutflow',
             'losPlane','losOutflow']

galNums = range(21,30)
a = 0.490
planeAzCut = 30
outflowAzCut = 40

for galNum in galNums:

    print('Galaxy = {0:d}'.format(galNum))

    loc = baseloc+subloc.format(galNum,a)
    fname = loc+cellFileName.format(galNum,a)
    rname = loc+rotFileName.format(a)
    lname = loc+linesName
    oname = outFileName.format(galNum,a)
    
    try:
        cells = pd.read_hdf(fname,'data')
        rotmat = pd.read_csv(rname,sep='\s+')
        lines = pd.read_csv(lname,sep='\s+',skiprows=2,
                            names=linesHeader)
    except IOError as e:
        print(galNum)
        continue

    # Only include the relavant columns
    cells = cells[cols]

    # Rotate the box to the galaxy frame
    print('\tRotate')
    cells = cells.apply(rotate,args=[rotmat],axis=1)
    
    # Convert to spherical coordinates
    print('\tConvert')
    cells = cells.apply(spherical,axis=1)
    
    # Add in Rvir
    cells['rvir'] = rotmat['Rvir']
    
    # Determine if the cell is in the plane or in outflows
    print('\tLabel')
    cells['cellPlane'] = cells['phi'].apply(lambda x: (90-planeAzCut<=x) & (x<=90.+planeAzCut))
    cells['cellOutflow'] = cells['phi'].apply(lambda x: (90-planeAzCut<=x) & (x<=90.+planeAzCut))

    #cells['cellPlane'] = cells['phi'].apply(cell_plane_cut,args=[planeAzCut])
    #cells['cellOutflow'] = cells['phi'].apply(cell_outflow_cut,args=[outflowAzCut])
    
    # Determine if the LOS is in the plane or in outflows
    #cells['losPlane'] = cells['LOS'].apply(los_plane_cut,args=[planeAzCut])
    #cells['losOutflow'] = cells['LOS'].apply(los_outflow_cut,args=[outflowAzCut])
    cells['losPhi'] = cells['LOS'].apply(lambda x: lines['phi'].iloc[int(x)-1])
    cells['losAz'] = cells['losPhi'].apply(azimuthal)
    cells['losPlane'] = cells['losAz'].apply(lambda x: x<=planeAzCut)
    cells['losOutflow'] = cells['losAz'].apply(lambda x: x>=outflowAzCut)
    
    print('\tSave')
    cells.to_hdf(oname,'data',mode='w')
    

    
    



