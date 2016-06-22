
import numpy as np
import pandas as pd


baseLoc = '/mnt/cluster/abs/cgm/'
baseGals = ['vela2a', 'vela2a']
galNums = [25,26,27,28]
expn = '0.490'
header = ['cellID', 'xbox', 'ybox','zbox','xgal','ygal','zgal','xsky','ysky','zsky']

incdeg = 90
inc = np.radians(incdeg)

# Sky to galaxy
stg_a11 = 1.
stg_a12 = 0.
stg_a13 = 0.

stg_a21 = 0.
stg_a22 = np.cos(inc)
stg_a23 = -1.0*np.sin(inc)

stg_a31 = 0.
stg_a32 = np.sin(inc)
stg_a33 = np.cos(inc)

# Galaxy to sky
gts_a11 = 1.
gts_a12 = 0.
gts_a13 = 0.

gts_a21 = 0.
gts_a22 = np.cos(inc)
gts_a23 = np.sin(inc)

gts_a31 = 0.
gts_a32 = -1.0*np.sin(inc)
gts_a33 = np.cos(inc)


for gal in baseGals:
    for galNum in galNums:

        loc = '{0:s}/{1:s}/vela{2:d}/a{3:s}/'.format(baseLoc,gal,galNum,expn)
        boxfile = loc+'{0:s}-{1:d}_GZa{2:s}.txt'.format(gal,galNum,expn)
        rotmatfile = loc+'rotmat_a{0:s}.txt'.format(expn)

        print boxfile

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

        coord = np.zeros((len(xbox),10))
        for i in range(len(xbox)):
           
            # Convert to galaxy frame
            xgal = a11*xbox[i] + a12*ybox[i] + a13*zbox[i]
            ygal = a21*xbox[i] + a22*ybox[i] + a23*zbox[i]
            zgal = a31*xbox[i] + a32*ybox[i] + a33*zbox[i]

            # Convert to sky frame
            xsky = gts_a11*xgal + gts_a12*ygal + gta_13*zgal 
            ysky = gts_a21*xgal + gts_a22*ygal + gta_23*zgal 
            zsky = gts_a31*xgal + gts_a32*ygal + gta_33*zgal

            coord[i,0] = i+1
            coord[i,1] = xbox[i]
            coord[i,2] = ybox[i]
            coord[i,3] = zbox[i]
            coord[i,4] = xgal
            coord[i,5] = ygal
            coord[i,6] = zgal
            coord[i,7] = xsky
            coord[i,8] = ysky
            coord[i,9] = zsky
        
        # Write to file
        outfile = '{0:s}/cell_i{1:d}.coords'.format(loc,incdeg)
        df = pd.DataFrame(coord, columns=header)
        df.to_hdf(outfile, 'data', mode='w')
     














