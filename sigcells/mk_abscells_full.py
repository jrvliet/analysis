
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
    form = '{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\t{5:s}\t{6:s}\t{7:s}\t{8:s}\n'
    for absf in absfiles:
        print '\tAbsf =',absf
        galID = absf.split('.')[0].split('/')[-1]
        
        cellid, expn = np.loadtxt(absf,skiprows=1, usecols=(2,13), 
                                    unpack=True)

        numcells = len(cellid)

        # Select 1000 random cells from this file
        for i in range(0,numcells):

            #ind = np.random.randint(0,numcells)
            ind = i
            a = '{0:.3f}'.format(expn[ind])
            cellnum = cellid[ind]

            # Get the cell properties
            
            gasfile = '{0:s}/{1:s}_GZa{2:s}.{3:s}.txt'.format(gasLoc,galID,a,ion)
            with open(gasfile, 'r') as f:

                for j,line in enumerate(f):
                    if j==cellnum+1:
                        l = line.split()
                        x = l[1]
                        y = l[2]
                        z = l[3]
                        vx = l[4]
                        vy = l[5]
                        vz = l[6]
                        nH = l[7]
                        temp = l[8]
                        snII = l[9]
                        fout.write(form.format(x,y,z,vx,vy,vx,nH,temp,snII))
                        break
    fout.close()
                                    
