#
# Plots the phase diagram for each ion/feedback method
# Each bin contains the total mass for that ion in the simulation
# Only cells that contribute significantly to absorption are considered
#
# To only plot the ISM, run as
#    python nT_mass.py ism
# To only plot the CGM, run as
#    python nT_mass.py cgm

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import linecache as lc


def getBin( value, bins ):
    index = -1
    for i in range(0, len(bins)-1 ):
        if value>bins[i] and value<=bins[i+1]:
            index = i

    return index

def printHisto( H, numbins ):

    f = open('histo.output', 'w')
    for i in range(0, numbins):
        for j in range(0, numbins):
            f.write(H[i][j] )
        f.write('\n')
    f.close()

def convertDensity( numDense, ion ):

    # Converts number density of an ion to mass density in 
    # solar masses per cubic parsec
    
    # Constants
    amu = 1.6605e-24  # in grams
    solar = 1.989e33  # in grams
    pc = 3.0856e18     # in cm

    atomicMass = getAtomicMass( ion )
    
    # Unlog the number density
#    numDense = pow(10, numDense)

    # Convert to solar mass
    density = numDense * atomicMass * amu / solar

    # Convert to cubic parsec
    density = density * pow(pc,3)

    return density


def getAtomicMass( ion ):
    mass = 0.0
    if ion=='HI':
        mass = 1.008
    elif ion=='MgII':
        mass = 24.305
    elif ion=='CIV':
        mass = 12.0107
    elif ion=='OVI':
        mass = 15.9994

    return mass

#print tempBins, len(tempBins)
'''
dumDense, dumTemp = [], []
for i in range(0,1000):
    dumDense.append( (maxDense - minDense) * np.random.random() + minDense )
    dumTemp.append( (maxTemp - minTemp) * np.random.random() + minTemp )
histos = []
'''
def inRegion( region, r):
    # Checks if the cell should be used
    if region=='ism':
        if r<0.1*rvir[i]:
            use = 1
        else:
            use = 0
    elif region=='cgm':
        if r>0.1*rvir[i]:
            use = 1
        else:
            use = 0
    elif region=='all':
        use = 1
    else:
        print 'Unknown region {0:s}'.format(region)
        sys.exit()
    return use

# Add the ability to split the CGM from the ISM
# Define the ISM as all cells with 0.1 Rvir of the galaxy
if len(sys.argv)>1:
    region = sys.argv[1]
else:
    region = 'all'



galID_list = ['D9o2', 'D9q', 'D9m4a']
feedback_list = ['dwSN', 'dwALL\_1', 'dwALL\_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']
rvir = [80.25497, 82.08840, 79.05807]

numbins = 50

dataLoc = '/home/matrix3/jrvander/sebass_gals/dwarfs/abscells/'
#dataLoc = '/home/jacob/matrix/sebass_gals/dwarfs/abscells/'

minDense = -8.0
maxDense = 1.0
minTemp = 2.0
maxTemp = 8.0

denseBins = np.linspace(minDense, maxDense, num=numbins+1, endpoint=True)
tempBins = np.linspace(minTemp, maxTemp, num=numbins+1, endpoint=True)
#print denseBins, len(denseBins)
#print tempBins, len(tempBins)

dumDense, dumTemp = [], []
for i in range(0,1000):
    dumDense.append( (maxDense - minDense) * np.random.random() + minDense )
    dumTemp.append( (maxTemp - minTemp) * np.random.random() + minTemp )
histos = []

fo = open('mass{0:s}.out'.format(region), 'w')
for i in range(0,len(galID_list)):
    galID = galID_list[i]
    galTotal = 0. 
    for ion in ion_list:
        
        # Initialize an empty 2D array
        H = np.zeros((numbins,numbins))
 
        # Open the data file
        filename = galID+'.'+ion+'.bulk_abscells_extend.dat'
        f = open(dataLoc+filename)
        f.readline()   # Read past header

        for line in f:
            
            l = line.split()
            r = float( l[6] )
            nH = float( l[7] )
            t = float( l[8] )
            cellsize = float( l[9] )  # in kpc
            mf = float( l[10] )
            alphaZ = float( l[12] )        
            cellID = int( l[2] )
            nion = float( l[14] )
            expn = float( l[13] )

            # Check if the cell should be used
            use = inRegion( region, r)
            # Get the index of this cell
            if use==1:
                denseInd = getBin(nH, denseBins)
                tempInd = getBin(t, tempBins)
                ionDense = convertDensity( nion, ion ) 
                cellsize = cellsize * 1000.    # Convert to pc from kpc
                cellMass = ionDense*math.pow(cellsize,3)
                H[denseInd][tempInd] += cellMass
            
        totalmass = np.sum(H)
        galTotal += totalmass
        fo.write('Total Mass in {0:s}, {1:s}: {2:.3e}\n'.format(galID, ion, totalmass))
        H = np.ma.masked_where(H==0, H)
        H = np.ma.log10( H )
#        H = np.log10(H)
        H = np.rot90(H)
        H = np.flipud(H)
        histos.append(H)

    fo.write('Total Mass in {0:s}: {1:.3e}\n\n'.format(galID, galTotal))

fo.close()


# Plot 
print 'Start Plotting'
labels = ['(a)', '(d)', '(g)', '(j)', 
          '(b)', '(e)', '(h)', '(k)', 
          '(c)', '(f)', '(i)', '(l)']

fig,((p11,p12,p13), (p21,p22,p23), (p31,p32,p33), (p41,p42,p43)) = plt.subplots(4,3,figsize=(10.2,10.8))
plot_list = [p11, p21, p31, p41, p12, p22, p32, p42, p13, p23, p33, p43]

for i in range(0,len(plot_list)):

    ax = plot_list[i]
    H = histos[i]
#    print '\n'
#    print np.amax(H)


    if i==0 or i==4 or i==8:
        # Ion is HI
        cbarmin = -5
        cbarmax = 5
        cbarLabel = '$\log$ ($M_{HI} / M_{\odot}$)'
    elif i==1 or i==5 or i==9:
        # Ion is MgII
        cbarmin = -5
        cbarmax = 2
        cbarLabel = '$\log$ ($M_{MgII} / M_{\odot}$)'
    elif i==2 or i==6 or i==10:
        # Ion is CIV
        cbarmin = -5
        cbarmax = 2
        cbarLabel = '$\log$ ($M_{CIV} / M_{\odot}$)'
    elif i==3 or i==7 or i==11:
        cbarmin = -5
        cbarmax = 2
        cbarLabel = '$\log$ ($M_{OVI} / M_{\odot}$)'
    ticks = np.arange(cbarmin, cbarmax+1)
    mesh = ax.pcolormesh(denseBins, tempBins, H, vmin=cbarmin, vmax=cbarmax)
    ax.set_xlim( [-8,1] )
    ax.set_ylim( [2,8] )
    ax.text( -1, 7, labels[i] )

#    if i==10:
#        for j in range(0,len(H[0])):
#            for k in range(0,len(H[0])):
#                n = denseBins[j]
#                t = tempBins[k]
#                if n>-6 and n<-3 and t>4 and t<5:
#                    print denseBins[j], tempBins[k], H[i][k]
#        printHisto(H, numbins)
     
    if i==0:
        ax.set_title('dwSN', size=12)
        ax.set_ylabel('HI \n $\log$ (T) [K] ', multialignment='center')
    elif i==1:
        ax.set_ylabel('MgII \n $\log$ (T) [K] ', multialignment='center')
    elif i==2:
        ax.set_ylabel('CIV \n $\log$ (T) [K] ', multialignment='center')
    elif i==3:
        ax.set_ylabel('OVI \n $\log$ (T) [K] ', multialignment='center')
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==4:
        ax.set_title('dwALL\_1', size=12)
    elif i==7:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    elif i==8:
        ax.set_title('dwALL\_8', size=12)
    elif i==11:
        ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    
    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(cbarLabel, rotation=270)
    cbar.set_ticks(ticks)

dropx = [p12, p22, p32, p13, p23, p33, p11, p21, p31]
plt.setp([a.get_xticklabels() for a in dropx],visible=False)

dropy = [p12, p22, p32, p13, p23, p33, p42, p43]
plt.setp([a.get_yticklabels() for a in dropy],visible=False)


plt.tight_layout()
s = dataLoc+'abscell_phase_master_mass{0:d}_{1:s}.pdf'.format(numbins,region)
plt.savefig(s, bbox_inches='tight')




