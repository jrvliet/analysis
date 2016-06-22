
import numpy as np
import pandas as pd


galNums = [25,26,27,28]
header = ['xbox', 'ybox','zbox','xgal','ygal','zgal','xsky','ysky','zsky']

inc = np.radians(90)

stg_a11 = 1.
stg_a12 = 0.
stg_a13 = 0.

stg_a21 = 0.
stg_a22 = np.cos(inc)
stg_a23 = -1.0*np.sin(inc)

stg_a31 = 0.
stg_a32 = np.sin(inc)
stg_a33 = np.cos(inc)



for galNum in galNums:

    boxfile = './vela{0:d}/a0.490/vela2b-{0:d}_GZa0.490.txt'.format(galNum)
    rotmatfile = './vela{0:d}/a0.490/rotmat_a0.490.txt'.format(galNum)

    with open(rotmatfile, 'r') as f:

        f.readline()
        l = f.readline().split()
        a11 = float(l[4])
        a12 = float(l[5])
        a13 = float(l[6])
        a21 = float(l[7])
        a22 = float(l[8])
        a23 = float(l[9])
        a31 = float(l[10])
        a32 = float(l[11])
        a33 = float(l[12])
    
    # Read in the box coordinates
    xbox, ybox, zbox = np.loadtxt(boxfile, skiprows=2, usecols=(1,2,3), unpack=True)

    # Convert to galaxy frame
    xgal, ygal, zgal = [], [], []
    for i in range(len(xbox)):
       
        x = a11*xbox[i] + a12*ybox[i] + a13*zbox[i]
        y = a21*xbox[i] + a22*ybox[i] + a23*zbox[i]
        z = a31*xbox[i] + a32*ybox[i] + a33*zbox[i]

        xgal.append(x)
        ygal.append(y)
        zgal.append(z)

        # Convert to sky frame
        





 
