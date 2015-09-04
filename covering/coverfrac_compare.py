#!/usr/bin/python
#
# Filename: coverfrac.py
# Version: 3 (formal)
#
# Author: Jacob Vander Vliet
# Date: 04/03/2014
#
# Compute the non-cumulative covering fraction
# Cut the impact space into annulli
# Add a loop over all EW cuts
# Add in errors
# Fixed the name of galaxies
#
# Plots the comparison across multiple feedback
# mechanisms for the same ion at the same redshift (set to z=0)
#
# Usage:
#  python coverfrac.py <ion> <ewcut>

def read_control_file(filename):

    f = open(filename)
    
    for i in range(0,5):
        line = f.readline()

    gasfile = line.split()[0]
    gasfile = gasfile.replace('"', '')
    galID = gasfile.split('_')[0]
    
    line = f.readline()
    rootname = line.split()[0]
    
    line = f.readline()
    aexpn = float(line.split()[0])
    
    for i in range(0,3):
        line = f.readline()
    summaryLoc = (line.split()[0]).replace('"','')

    f.close()

    return gasfile, galID, rootname, aexpn, summaryLoc

def read_summary(galID, aexpn, summaryLoc):

    import sys

    aexpn = float(aexpn)

    # Location of summaries
    f = open(summaryLoc + galID + '.dat')

    f.readline()
    f.readline()

    found = 0

    for line in f:
#    while(found==0):

        l = line.split()
        a = float(l[0])

        if a==aexpn:
            
            found = 1

            # Read in data
            mvir = float(l[2])
            rvir = float(l[3])
            a11 = float(l[4])
            a12 = float(l[5])
            a13 = float(l[6])
            a21 = float(l[7])
            a22 = float(l[8])
            a23 = float(l[9])
            a31 = float(l[10])
            a32 = float(l[11])
            a33 = float(l[12])
            vpec_x = float(l[13])
            vpec_y = float(l[14])
            vpec_z = float(l[15])
            
    f.close()
    if found==1:
        return mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z
    else:
        print 'ERROR: Could not find the expansion factor {0} in file: \n\t'.format(aexpn) + summaryLoc
        sys.exit()


# Function for determining vertical error bars
# Uses incomplete Beta functions
def funcup(x,n1,n2,cl):
    a = n1 + 1.0
    b = n2
    return cl - sc.betainc(a,b,x)

def funcdn(x,n1,n2,cl):
    a = n2 + 1.0
    b = n1
    return cl - sc.betainc(a,b,x)


import numpy as np
import scipy.special as sc
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import math
import matplotlib
import os
import subprocess as sp
import os.path as op


font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

#matplotlib.rc('font',**font)

plot_cos = 0

galID_list = ['dwRP_1', 'dwRP_8', 'dwSN']
galID_label = ['dwALL_1', 'dwALL_8', 'dwSN']
galID_file = ['q', 'm4a', 'o2']
ion = sys.argv[1]
expn_list = ['1.001', '1.000', '1.001']
mark=['v','o','x']
colors=['fuchsia', 'blue', 'black']
line=['-','.']

# Define the EW cuts to use
ewcut = float(sys.argv[2])

# Define the radial bins to use
nbins = 15
binsize = 1.5/nbins

i=0

################################################################################################################################


for i in range(0,len(galID_list)):

    # Move into that directory
    os.chdir('./D9'+galID_file[i]+'_outputs/a'+expn_list[i]+'/'+ion+'/')

    sp.call('pwd',shell=True)
    fo1 = open('covering1.out','w')
    fo2 = open('covering2.out','w')

    # Read in galaxy properties
    controlfile = 'gal_props_'+ion+'.dat'
    gasfile, galID, rootname, expn, summaryLoc = read_control_file(controlfile)
    galID = 'D9'+galID_file[i]

    mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z = read_summary(galID, expn, summaryLoc)

    # Determine the bin edges
    Dmax = []
    Dmin = [0.0]
    for j in range(0,nbins):
        Dmax.append( (j+1)*binsize )
        Dmin.append( (j+1)*binsize )
    del Dmin[-1]

    # COS-Dwarf bins
    if plot_cos==1:
        Dmin = [0.0, 0.25, 0.55]
        Dmax = [0.25, 0.55, 1.0]
    

    los=[]
    CL = 0.8413
    tol = 1e-4
    pmin = 0.0
    pmax = 1.0
    iter = 1000

    # Open file with list of all LOS
    f = open('lines.info','r')
    f.readline()
    f.readline()
    for line in f:
        los.append(line.split())
    f.close()

    index = 0
#    sp.call('pwd', shell=True)
    allfile = galID+'.'+ion+'.ALL.sysabs'
#    print ''
#    print allfile
    
    impact=[]
    abs=[]
    absimpact=[]
    covering=[]
    imp=[]
    horizerrneg=[]
    horizerrpos=[]
    verterrneg=[]
    verterrpos=[]


    if op.isfile(allfile):
        notthere = 0
        f = open(allfile,'r')
        header = (f.readline()).split()
        a = header[0]
        z = header[1]
        for line in f:
            abs.append(line.split())
        f.close()

    else:
        print 'No ALL file for '+galID+', '+ion
        notthere = 1
#        os.chdir('../../..')
#        break


                    
    # Loop over the different impact parameters
    for j in range(0,len(Dmax)):
        maxrad = Dmax[j]
        minrad = Dmin[j]
        hit = 0.0
        total = 0.0
        xaxispoint = 0.0
        zerobinabs = []
        zerobin = []
        
        if notthere == 0:
            # Loop over all lines from the .sysabs file to get the
            # number of lines with significant absorption (SL > 3)
            # with an impact parameter between minrad and maxrad
            for k in range(0,len(abs)):
                absimpact.append(float(abs[k][1]) / rvir)
                width = float(abs[k][5])
                if absimpact[k]>minrad and absimpact[k]<maxrad and width>ewcut:
                    hit += 1
                    xaxispoint += absimpact[k]
                    if minrad==0.0:
                        zerobinabs.append( abs[k][0] )
            
            
            # Loop over all lines from impactparam.dat to get the
            # total number of sightlines with an impact parameter
            # between minrad and maxrad
            for k in range(0,len(los)):
                impact.append(float(los[k][1]) / rvir)
                if impact[k] > minrad and impact[k] < maxrad:
                    total += 1
                    if minrad==0.0:
                        continue
#                        print k, len(los), abs[k][0]
#                        zerobin.append( abs[k][0] )
                    
            fraction = hit/total
            covering.append(fraction)

#            if len(zerobinabs) !=  len(zerobin) and j==0:
#                print set(zerobinabs)^set(zerobin)
    
            # Determine the x location of the point. If hit>0, then the location is the
            # average impact parameter in that bin. If hit=0, then the location is the
            # midpoint of the bin.
            if hit > 0.0:
                imp.append(xaxispoint/hit)  # location of point on x-axis
            else:
                imp.append((maxrad-minrad)/2.0 + minrad) 
            
            # Determine error bars
            #   Horizontal error bars are width of bin
            horizerrneg.append(imp[j]-minrad)
            horizerrpos.append(maxrad-imp[j])
        
            #   Vertical error bars are found using the incomplete beta function
            top = opt.brentq(funcup,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)
            bot = 1.0 - opt.brentq(funcdn,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)

            if fraction == 1.0:
                top = 1.0
            elif fraction == 0.0:
                bot = 0.0

            verterrpos.append(top-fraction)
            verterrneg.append(fraction-bot)

        elif notthere==1:
            imp.append((maxrad-minrad)/2.0 + minrad)
            covering.append(0.0)
            horizerrneg.append(imp[j]-minrad)
            horizerrpos.append(maxrad-imp[j])
            verterrpos.append(0.03)
            verterrneg.append(0.00000)            

        z = math.fabs(1.0/float(expn) - 1.0)

    plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos],label=galID_label[i],marker=mark[i],c=colors[i],linestyle='none')
        
    # Write to output file
    if j == 0:
        for k in range(0,len(Dmax)):
            fo1.write('{0}\t{1}\t{2}\t{3}\n'.format(ewcut,Dmax[k],imp[k],covering[k]))
    elif j == 1:
        for k in range(0,len(Dmax)):
            fo2.write('{0}\t{1}\t{2}\t{3}\n'.format(ewcut,Dmax[k],imp[k],covering[k]))

    index += 1
    i += 1

    # Move up a directory
    os.chdir('../../..')


#sp.call('pwd', shell=True)

# If ion is CIV, overplot COS-Dwarfs Data
if ion == 'CIV' and plot_cos==1:
    imp = [0.15, 0.375, 0.750]
    covering = [0.875, 0.475, 0.0]
    verterrneg = [0.125, 0.125, 0.0]
    verterrpos = [0.125, 0.125, 0.05]

    horizerrneg = [0.15, 0.125, 0.20]
    horizerrpos = [0.10, 0.175, 0.25]

    print len(imp), len(covering), len(horizerrneg), len(horizerrpos), len(verterrneg), len(verterrpos)

    plt.errorbar(imp, covering, xerr=[horizerrneg,horizerrpos], yerr=[verterrneg,verterrpos], label='COS-Dwarfs', marker='D', c='red', linestyle='none')


#plt.title('z = 0; '+ion+ '; EW Cut: {0:0.1f} $\AA$'.format(ewcut) )
#plt.title(ion+ '; EW Cut: {0:0.1f} m$\AA$'.format(ewcut*1000) )
#plt.text(0.78, 0.82, '$W \,\\geq\,0.1\\AA$', fontsize=15)

plt.xlabel('$\\rho$ / R$_{vir}$')
ylab = 'C$_f$ ['+ion+']'
plt.ylabel(ylab)
plt.ylim([-0.1,1.1])
plt.xlim([min(Dmin),1.0])
plt.legend(numpoints=1,loc=1,prop={'size':15}, frameon=False)
plt.gca().get_frame().set_linewidth(2) 



plt.tight_layout()
z = math.fabs(1.0/float(expn) - 1.0)
red = '{0:.2f}'.format(z)
plt.draw()
plt.savefig(ion+'_{0:0.1f}mA_covercompare.eps'.format(ewcut*1000))





'''

################################################################################################################################

# Move up a directory
#print os.getcwd()
os.chdir('../../../D9m4a_outputs/a1.000/'+ion+'/')
#print os.getcwd()

# Read in galaxy properties
controlfile = 'gal_props_'+ion+'.dat'
gasfile, galID, rootname, expn, summaryLoc = read_control_file(controlfile)
mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z = read_summary(galID, expn, summaryLoc)

# Determine the bin edges
Dmax = []
Dmin = [0.0]
for j in range(0,nbins):
    Dmax.append( (j+1)*binsize )
    Dmin.append( (j+1)*binsize )
del Dmin[-1]

los=[]
CL = 0.8413
tol = 1e-4
pmin = 0.0
pmax = 1.0
iter = 1000

# Open file with list of all LOS
f = open('lines.info','r')
f.readline()
f.readline()
for line in f:
    los.append(line.split())

index = 0

allfile = galID+'.'+ion+'.ALL.sysabs'
f = open(allfile,'r')

header = (f.readline()).split()

a = header[0]
z = header[1]

impact=[]
abs=[]
absimpact=[]
covering=[]
imp=[]
horizerrneg=[]
horizerrpos=[]
verterrneg=[]
verterrpos=[]
for line in f:
    abs.append(line.split())
f.close()
                    
# Loop over the different impact parameters
for j in range(0,len(Dmax)):
    maxrad = Dmax[j]
    minrad = Dmin[j]
    hit = 0.0
    total = 0.0
    xaxispoint = 0.0

    zerobinabs = []
    zerobin = []

    # Loop over all lines from the .sysabs file to get the
    # number of lines with significant absorption (SL > 3)
    # with an impact parameter between minrad and maxrad
    for k in range(0,len(abs)):
        absimpact.append(float(abs[k][1]) / rvir)
        width = float(abs[k][5])
        if absimpact[k]>minrad and absimpact[k]<maxrad and width>ewcut:
            hit += 1
            xaxispoint += absimpact[k]
            if minrad==0.0:
                zerobinabs.append( abs[k][0] )

    # Loop over all lines from impactparam.dat to get the
    # total number of sightlines with an impact parameter
    # between minrad and maxrad
    for k in range(0,len(los)):
        impact.append(float(los[k][1]) / rvir)
        if impact[k] > minrad and impact[k] < maxrad:
            total += 1
#            if minrad==0.0:
#                zerobin.append( abs[k][0] )   


    fraction = hit/total
            
    covering.append(fraction)
        
    if len(zerobinabs) !=  len(zerobin) and j==0:
        print set(zerobinabs)^set(zerobin)

            
    # Determine the x location of the point. If hit>0, then the location is the
    # average impact parameter in that bin. If hit=0, then the location is the
    # midpoint of the bin.
    if hit > 0.0:
        imp.append(xaxispoint/hit)  # location of point on x-axis
    else:
        imp.append((maxrad-minrad)/2.0 + minrad) 
            
    # Determine error bars
    #   Horizontal error bars are width of bin
    horizerrneg.append(imp[j]-minrad)
    horizerrpos.append(maxrad-imp[j])
        
    #   Vertical error bars are found using the incomplete beta function
    top = opt.brentq(funcup,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)
    bot = 1.0 - opt.brentq(funcdn,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)

    if fraction == 1.0:
        top = 1.0
    if fraction == 0.0:
        bot = 0.0

    verterrpos.append(top-fraction)
    verterrneg.append(fraction-bot)

z = math.fabs(1.0/float(expn) - 1.0)

plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos],label=galID_list[i],marker=mark[i],c=colors[i],linestyle='none')
        
# Write to output file
if j == 0:
    for k in range(0,len(Dmax)):
        fo1.write('{0}\t{1}\t{2}\t{3}\n'.format(ewcut,Dmax[k],imp[k],covering[k]))
elif j == 1:
    for k in range(0,len(Dmax)):
        fo2.write('{0}\t{1}\t{2}\t{3}\n'.format(ewcut,Dmax[k],imp[k],covering[k]))

index += 1
i += 1


################################################################################################################################

# Move to last file only is ion is not OVI or CIV
# This is becuase the SN only model fails to generate
# any OVI or CIV

if ion!='OVI' and ion!='CIV':

    os.chdir('../../../preliminary/dwarf9o/'+ion+'/aas/')

    
    # Read in galaxy properties
    controlfile = 'gal_props_'+ion+'.dat'
    gasfile, galID, rootname, expn, summaryLoc = read_control_file(controlfile)
    mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z = read_summary(galID, expn, summaryLoc)
    
    # Determine the bin edges
    Dmax = []
    Dmin = [0.0]
    for j in range(0,nbins):
        Dmax.append( (j+1)*binsize )
        Dmin.append( (j+1)*binsize )
    del Dmin[-1]

    los=[]
    CL = 0.8413
    tol = 1e-4
    pmin = 0.0
    pmax = 1.0
    iter = 1000
    
    # Open file with list of all LOS
    f = open('lines.info','r')
    f.readline()
    f.readline()
    for line in f:
        los.append(line.split())

    index = 0

    allfile = galID+'.'+ion+'.ALL.sysabs'
    f = open(allfile,'r')
    
    header = (f.readline()).split()
    
    a = header[0]
    z = header[1]
    
    impact=[]
    abs=[]
    absimpact=[]
    covering=[]
    imp=[]
    horizerrneg=[]
    horizerrpos=[]
    verterrneg=[]
    verterrpos=[]
    for line in f:
        abs.append(line.split())
    f.close()
                    
    # Loop over the different impact parameters
    for j in range(0,len(Dmax)):
        maxrad = Dmax[j]
        minrad = Dmin[j]
        hit = 0.0
        total = 0.0
        xaxispoint = 0.0

        zerobinabs = []
        zerobin = []

        # Loop over all lines from the .sysabs file to get the
        # number of lines with significant absorption (SL > 3)
        # with an impact parameter between minrad and maxrad
        for k in range(0,len(abs)):
            absimpact.append(float(abs[k][1]) / rvir)
            width = float(abs[k][5])
            if absimpact[k]>minrad and absimpact[k]<maxrad and width>ewcut:
                hit += 1
                xaxispoint += absimpact[k]
                if minrad==0.0:
                    zerobinabs.append( abs[k][0] )   

        # Loop over all lines from impactparam.dat to get the
        # total number of sightlines with an impact parameter
        # between minrad and maxrad
        for k in range(0,len(los)):
            impact.append(float(los[k][1]) / rvir)
            if impact[k] > minrad and impact[k] < maxrad:
                total += 1
                if minrad==0.0:
                    zerobin.append( abs[k][0] )   


        fraction = hit/total
            
        covering.append(fraction)
                 

        if len(zerobinabs) !=  len(zerobin) and j==0:
            print set(zerobinabs)^set(zerobin)

   
        # Determine the x location of the point. If hit>0, then the location is the
        # average impact parameter in that bin. If hit=0, then the location is the
        # midpoint of the bin.
        if hit > 0.0:
            imp.append(xaxispoint/hit)  # location of point on x-axis
        else:
            imp.append((maxrad-minrad)/2.0 + minrad) 
            
        # Determine error bars
        #   Horizontal error bars are width of bin
        horizerrneg.append(imp[j]-minrad)
        horizerrpos.append(maxrad-imp[j])
        
        #   Vertical error bars are found using the incomplete beta function
        top = opt.brentq(funcup,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)
        bot = 1.0 - opt.brentq(funcdn,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)

        if fraction == 1.0:
            top = 1.0
        if fraction == 0.0:
            bot = 0.0

        verterrpos.append(top-fraction)
        verterrneg.append(fraction-bot)

    z = math.fabs(1.0/float(expn) - 1.0)

    plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos],label=galID_list[i],marker=mark[i],c=colors[i],linestyle='none')
        
    # Write to output file
    if j == 0:
        for k in range(0,len(Dmax)):
            fo1.write('{0}\t{1}\t{2}\t{3}\n'.format(ewcut,Dmax[k],imp[k],covering[k]))
    elif j == 1:
        for k in range(0,len(Dmax)):
            fo2.write('{0}\t{1}\t{2}\t{3}\n'.format(ewcut,Dmax[k],imp[k],covering[k]))

    index += 1


    # Move back to original directory
    os.chdir('../../../../')

else:
    os.chdir('../../../')
    imp = []
    verterrpos = []
    verterrneg = []
    covering = []
    horizerrpos = []
    horizerrneg = []
    for j in range(0,len(Dmax)):
        maxrad = Dmax[j]
        minrad = Dmin[j]
        imp.append((maxrad-minrad)/2.0 + minrad) 
        covering.append(0.0)
        verterrpos.append(0.0)
        verterrneg.append(0.0)
        horizerrneg.append(imp[j]-minrad)
        horizerrpos.append(maxrad-imp[j])
    
    plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos],label=galID_list[i],marker=mark[i],c=colors[i],linestyle='none')




#plt.title('z = 0; '+ion+ '; EW Cut: {0:0.1f} $\AA$'.format(ewcut) )
plt.title(ion+ '; EW Cut: {0:0.1f} m$\AA$'.format(ewcut*1000) )
plt.xlabel('Impact Parameter (Rvir)')
plt.ylabel('Covering Fraction')
plt.ylim([-0.1,1.1])
plt.xlim([min(Dmin),max(Dmax)])
plt.legend(numpoints=1,loc=1,prop={'size':10}, frameon=False)
plt.gca().get_frame().set_linewidth(2) 



plt.tight_layout()
z = math.fabs(1.0/float(expn) - 1.0)
red = '{0:.2f}'.format(z)
plt.draw()
plt.savefig(ion+'_{0:0.1f}mA_covercompare.pdf'.format(ewcut*1000))
'''
