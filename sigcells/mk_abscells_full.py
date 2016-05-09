
'''
Geneartes subsets of the <ion>_vlos.dat files small enough
for heirarchical clustering testing

Output files are names <ion>_vlos.subset and are no more than
1000 lines long

'''
import sys
import numpy as np
from glob import glob

ions = ['HI', 'MgII', 'CIV', 'OVI']
subsetSize = 1000

absLoc = '/home/jacob/research/dwarfs/abscells/'
absLoc = '/home/hyades/jrvander/dwarfs/abscells/'
gasLoc = absLoc.replace('abscells', 'gasfiles')


for ion in ions:

    print ion
    # Get the list of files of absorbing cells
    absfiles = glob('{0:s}*.{1:s}.bulk_abscells.dat'.format(
                    absLoc, ion))

    # Open the output file
    outfile = '{0:s}_abscells.full'.format(ion)
    fout = open(outfile, 'w')
    fout.write('x\t\ty\t\tz\t\tvx\t\tvy\t\tvz\t\tnH\t\tTemp\t\tsnII\n')
    form = '{0:.4e}\t{1:.4e}\t{2:.4e}\t{3:.4e}\t{4:.4e}\t{5:.4e}\t{6:.4e}\t{7:.4e}\t{8:.4e}\n'
    for absf in absfiles:
        print '\tAbsf =',absf
        galID = absf.split('.')[0].split('/')[-1]
        
        cellid, expn = np.loadtxt(absf,skiprows=1, usecols=(2,13), 
                                    unpack=True)

        numcells = len(cellid)

        prevA = '0.0'
        # Select 1000 random cells from this file
        for i in range(0,numcells):

            #ind = np.random.randint(0,numcells)
            ind = i
            a = '{0:.3f}'.format(expn[ind])
            cellnum = cellid[ind]

            # If a is a different one from before, read in a new gas box
            if a!=prevA:
                print 'New expansion factor: {0:s} -> {1:s}'.format(prevA, a)
                gasfile = '{0:s}/{1:s}_GZa{2:s}.{3:s}.txt'.format(gasLoc,galID,a,ion)
                gas = np.loadtxt(gasfile,skiprows=2)
                prevA = a

            index = int(cellnum-1)
            x = gas[index,1]
            y = gas[index,2]
            z = gas[index,3]
            vx = gas[index,4]
            vy = gas[index,5]
            vz = gas[index,6]
            nH = gas[index,7]
            temp = gas[index,8]
            snII = gas[index,9]
            
            fout.write(form.format(x,y,z,vx,vy,vx,nH,temp,snII))


    fout.close()
