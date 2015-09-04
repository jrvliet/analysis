#!/usr/bin/python
#
# Filename: plotWvD_master.py
# Date: 04/09/14
#
# Usage:
#  python plotWvD_master.py <ewcut in Angstroms> <bulk?>
#
# Plots the time evolution of a single ion
# Compares the different feedback perscriptions
# Bulk flag, if =1, will plot all four ions on one plot


import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv)!=3:
    bulk = 0
else:
    bulk = sys.argv[2]


galID_list = ['D9o2', 'D9q', 'D9m4a']
labels = ['dwSN', 'dwALL\_1', 'dwALL\_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']
mark=['v','o','x']
colors=['fuchsia', 'blue', 'black']
line=['-','.']

ewcut = float(sys.argv[1])
record = 1

# Define the radial bins to use
nbins = 15
binsize = 1.5/nbins

Dmax = []
Dmin = [0.0]
for j in range(0,nbins):
    Dmax.append( (j+1)*binsize )
    Dmin.append( (j+1)*binsize )
del Dmin[-1]
    
for ion in ion_list:
    outputs = []
    for i in range(0,len(galID_list)):
        galID = galID_list[i]

        # Read in the master ALL file
        master_file = galID+'.'+ion+'.ALL.sysabs.large.master'
        data = np.loadtxt(master_file,skiprows=1)
        
        # Loop over the different impact parameter bins
        imp = []
        horizerrneg=[]
        horizerrpos=[]
        verterrneg=[]
        verterrpos=[]
        ewmean = []
        
        for j in range(0,len(Dmax)):
            minrad = Dmin[j]
            maxrad = Dmax[j]
            
            hit = 0.0
            xaxispoint = 0.0
            ewtotal = []

            # Loop over all the lines in the .sysabs file to get
            # the number of lines with impact parameters in this 
            # bin
            for k in range(0,len(data[:,0])):
                impact = data[k,22]
                ew = data[k,5]
                ew_err = data[k,6]
                if impact<maxrad and impact>minrad and ew_err!=-1.0 and ew>ewcut:
                    hit += 1
                    xaxispoint += impact
                    ewtotal.append(ew)

            # Take the mean of ewtotal
            ewmean.append(np.mean(ewtotal))
    

            # Determine the x location of the point. If hit>0, then the location is
            # the average impact parameter in that bin. If hit=0, then the location
            # is the midpoint of the bin
            if hit > 0.0:
                imp.append(xaxispoint/hit)  # location of point on x-axis
            else:
                imp.append((maxrad-minrad)/2.0 + minrad)

            # Determine error bars
            #   Horizontal error bars are width of bin
            horizerrneg.append(imp[j]-minrad)
            horizerrpos.append(maxrad-imp[j])
           
            #   Vertical error bars are found using the standard deviation
            standdev = np.std(ewtotal)
            if standdev<=np.mean(ewtotal):
                verterrpos.append(standdev)
                verterrneg.append(standdev)
            else:
                verterrpos.append(standdev)
                verterrneg.append(np.mean(ewtotal)-1e-5)

#    imp = [k/rvir for k in imp]
#    horizerrneg = [k/rvir for k in horizerrneg]
#    horizerrpos = [k/rvir for k in horizerrpos]

        # Write results to file
        if record==1:
            fout = open(galID+'.'+ion+'.'+'cut_{0:d}mA_ewImp.out'.format(int(ewcut*1000)), 'w')
            fout.write('d/Rvir \t <EW> \t +err \t -err \t flag\n')
            for k in range(0,len(imp)):
                if verterrpos[k]>ewmean[k]:
                    flag=1
                else:
                    flag=0
                fout.write('{0:.3f} \t {1:.3f} \t {2:.3f} \t {3:.3f} \t {4:d}\n'.format(imp[k], ewmean[k], verterrpos[k], verterrneg[k], flag))
            fout.close()
        
    
        if bulk==0:
            plt.errorbar(imp, ewmean, xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos], c=colors[i], marker=mark[i], label=labels[i], linestyle='none')
        else:
            subplotnum = 221+ion_list.index(ion)
            plt.subplot(subplotnum)
            plt.errorbar(imp, ewmean, xerr=[horizerrneg,horizerrpos],yerr=[verterrneg,verterrpos], c=colors[i], marker=mark[i], label=labels[i], linestyle='none')


    if bulk==0:
        plt.xlabel('D / R$_{vir}$')
        ylab = 'EW$_{'+ion+'}$ [$\AA$]'
        plt.ylabel(ylab)
        plt.yscale('log')
        ymin, ymax = plt.ylim()
        plt.xlim([-0.05,1.5])
        if ion=='HI':
            plt.ylim([1e-2, 20.0])
        else:
            plt.ylim([5e-4, 2.0])
        plt.legend(numpoints=1,loc=1,prop={'size':15}, frameon=False)
        
        s = 'master_'+ion+'_ewspace_cut_{0:d}mA.eps'.format(int(ewcut*1000))
        print s
        
        plt.savefig(s)
        
        s = 'master_'+ion+'_ewspace_cut_{0:d}mA.pdf'.format(int(ewcut*1000))
        plt.savefig(s)
        
        s = 'master_'+ion+'_ewspace_cut_{0:d}mA.jpg'.format(int(ewcut*1000))
        plt.savefig(s)


        plt.cla()
        plt.clf()


    else:
        subplotnum = 221+ion_list.index(ion)
        plt.subplot(subplotnum)
        plt.xlabel('D / R$_{vir}$')
        ylab = 'EW$_{'+ion+'}$ [$\AA$]'
        plt.ylabel(ylab)
        plt.yscale('log')
        ymin, ymax = plt.ylim()
        plt.xlim([-0.05,1.5])
        if ion=='HI':
            plt.ylim([1e-2, 20.0])
        else:
            plt.ylim([5e-4, 2.0])
        plt.legend(numpoints=1,loc=1,prop={'size':8}, frameon=False)
        
        


if bulk!=0:
    plt.subplots_adjust(wspace=0.3)
    s = 'master_ewspace_cut_{0:d}mA_bulk.eps'.format(int(ewcut*1000))
    print s
    
    plt.savefig(s, bbox_inches='tight')
    
