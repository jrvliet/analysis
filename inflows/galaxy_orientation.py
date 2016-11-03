
'''
Determines the angle between the galaxy's orientation
as defined by the rotation matix and the vector along 
the inflow at z=1
'''

from __future__ import print_function
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


def unit_vector(vector):
    return vector / la.norm(vector)

def angle_between(v1,v2):
    '''
    Returns angle between vectors in radians
    '''
    v1u = unit_vector(v1)
    v2u = unit_vector(v2)
    dot = v1u[0]*v2u[0] + v1u[1]*v2u[1] + v1u[2]*v2u[2]  
    return np.arccos(dot)
    #return np.arccos(np.clip(np.dot(v1u.T,v2u),-1,-1))

# Read in inflow matrix
loc = '/home/jacob/research/code/analysis/inflows/'
filename = 'inflow_rot_mat.dat'
aInflow = np.matrix(np.loadtxt(loc+filename))
aInflowInv = la.inv(aInflow)

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
rotname = 'a{0:.3f}/rotmat_a{0:.3f}.txt'

expns = np.arange(0.200,0.500,0.01)
header = 'Expn\tGalz0\tGalz1\tGalz3\tGIDegrees\tBIDegrees\tBGDegrees'
numSnaps = len(expns)
numHeads = len(header.split('\t'))
print(numSnaps,numHeads)
angles = np.zeros((numSnaps,numHeads))

for i,expn in enumerate(expns):

    rname = dataloc+rotname.format(expn)

    with open(rname) as f:
        f.readline()
        l = f.readline().split()
        a = np.array(map(float,l[4:-3])).reshape(3,3)
        
    a = np.matrix(a)
    aInv = la.inv(a)
    
    # Get the zaxis of the the box coordinate system
    boxz = np.matrix([0,0,1])

    # Get the zaxis of the the galaxy coordinate system
    galz = aInv*boxz.T

    # Get the zaxis of the inflow coordinate system
    inflowz = aInflowInv*boxz.T
    

    # Get the angle between the two vectors
    giAngle = np.degrees(angle_between(galz,inflowz))
    gbAngle = np.degrees(angle_between(galz,boxz.T))
    biAngle = np.degrees(angle_between(boxz.T,inflowz))

    angles[i,0] = expn
    angles[i,1] = galz[0]
    angles[i,2] = galz[1]
    angles[i,3] = galz[2]
    angles[i,4] = giAngle
    angles[i,5] = biAngle       
    angles[i,6] = gbAngle
    

outFile = 'inflow_galaxy_angles.dat'
header = header.replace('\t','\t\t')
np.savetxt(outFile,angles,fmt='%.6f',delimiter='\t',header=header)


fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(angles[:,0],angles[:,4],label='Gal-Inflow')
ax.plot(angles[:,0],angles[:,5],label='Box-Inflow')
ax.plot(angles[:,0],angles[:,6],label='Gal-Box')
ax.legend()
ax.set_xlabel('Expn')
ax.set_ylabel('Angle [deg]')
fig.savefig('inflow_galaxy_angle.png',bbox_inches='tight',dpi=300)


