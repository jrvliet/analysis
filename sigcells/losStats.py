
import numpy as np
import matplotlib.pyplot as plt


ions = []

with open('mockspec.config', 'r') as f:
    for i in range(23):
        f.readline()

    for line in f:
        ions.append(line.split()[0])
print ions

with open('galaxy.props', 'r') as f:
    galID = f.readline().split()[0]
    expn = f.readline().split()[0]
    redshift = f.readline().split()[0]
    mvir = f.readline().split()[0]
    rvir = f.readline().split()[0]
    inc = f.readline().split()[0]

with open('los_stats.txt', 'w') as f:
    
    s = ''
    for ion in ions:
        s += '\t\t\t{0:s}'.format(ion)
    f.write(s+'\n')
    


    for i in range(1000):
        losnum = str(i).zfill(4)
        origNum = []
        finNum = []
    
        s = '{0:d}\t'.format(i)

        for ion in ions:
            linesfile = '{0:s}/{1:s}.{0:s}.los{2:s}.lines'.format(ion, galID, losnum)
            final = '{0:s}.final'.format(linesfile)
    
            # Get the number of cells
            num = 0
            with open(linesfile, 'r') as flines:
                flines.readline()
                for line in flines:
                    num += 1        

            finalnum = 0
            with open(final, 'r') as flines:
                flines.readline()
                for line in flines:
                    finalnum += 1        
            
            percent = float(finalNum)/float(num)*100
            s += '{0:d}\t{1:d}\t{2:.2f}'.format(num,finalnum,percent)
        
        # Write to file
        f.write(s+'\n')


        
    





