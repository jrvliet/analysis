
import os.path
import numpy as np
import matplotlib.pyplot as plt


ions = []

with open('mockspec.config', 'r') as f:
    for i in range(24):
        f.readline()

    for line in f:
        ions.append(line.split()[0])

with open('galaxy.props', 'r') as f:
    galID = f.readline().split()[1]
    expn = f.readline().split()[1]
    redshift = f.readline().split()[1]
    mvir = f.readline().split()[1]
    rvir = f.readline().split()[1]
    inc = f.readline().split()[1]


with open('los_stats.txt', 'w') as f:
    
    s = '\t'
    for ion in ions:
        s += '{0:s}\t\t\t\t'.format(ion)
    f.write(s+'\n')
    

    origNum = np.zeros(len(ions))
    finNum = np.zeros(len(ions))
    count = np.zeros(len(ions))
    for i in range(1,1000):
        losnum = str(i).zfill(4)
    
        s = '{0:d}\t'.format(i)
        origSum = 0
        finalSum = 0
        ionnum = 0
        for ion in ions:
            linesfile = './{0:s}/{1:s}.{0:s}.los{2:s}.lines'.format(ion, galID, losnum)
            final = '{0:s}.final'.format(linesfile)
    
            # Get the number of cells
            num = 0
            with open(linesfile, 'r') as flines:
                flines.readline()
                for line in flines:
                    num += 1        

            finalnum = 0
            try:
                with open(final, 'r') as flines:
                    flines.readline()
                    for line in flines:
                        finalnum += 1        
            except IOError:
                finalnum = 0
            
            if num > 0:
                sysabs = './{0:s}/{1:s}.{0:s}.los{2:s}.sysabs'.format(ion, galID, losnum)
                if os.path.isfile(sysabs):
                    count[ionnum] += 1.0
                    origNum[ionnum] += num
                    finNum[ionnum] += finalnum            
            
            ionnum += 1
            try:
                percent = float(finalnum)/float(num)*100
            except ZeroDivisionError:
#                print 'Division error for los{0:d} in {1:s}: num={2:d}, finalnum={3:d}'.format(i,ion,num,finalnum)
                percent = 0
            s += '{0:d}\t{1:d}\t{2:.2f}\t|\t'.format(num,finalnum,percent)
        
        # Write to file
        f.write(s+'\n')

    print origNum
    print finNum

    s = ''
    for i in range(0,len(ions)*15):
        s += '-'
    f.write(s+'\n')

    s = 'Mean \t'
    for i in range(0,len(ions)):
        aveNum = float(origNum[i]) / count[i]
        aveFin =  float(finNum[i]) / count[i]
        try:
            aveper = aveFin / aveNum * 100
        except ZeroDivisionError:
            print ions[i]
            aveper = 0
        s += '{0:.1f}\t{1:.1f}\t{2:.2f}\t|\t'.format(aveNum,aveFin,aveper)
    f.write(s+'\n')


    





