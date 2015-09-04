#!/usr/bin/python
#
# Filename: ewImpBin_histogram.py
# Date: 16/12/14
# 
# Usage: 
#  python ewImpBin_histogram.py <ewcut in Angstrom> <bulk?>
#
# Plots a histogram of the EW vs impact parameter/Rvir
# Applys an equivilant width cut provided in command line
# Bulk flag, if =1, will plot all four ions on one plot

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def fmt2(x, pos):
    a = '{0:.0f}%'.format(x*100)
    return a


# Read in command line args
ewcut = float(sys.argv[1])
record = 1

if len(sys.argv)!=3:
    bulk = 0
else:
    bulk = sys.argv[2]


# Running Parameters
normalize = 0   # Normalize the histograms
do_logEW = 1    # determine if the EW should be logged
do_logCount = 0 # determine if the final bins should be logged

# Base of file names to use
galID_list = ['D9o2', 'D9q', 'D9m4a']
labels = ['dwSN', 'dwALL\_1', 'dwALL\_8']
ion_list = ['HI', 'MgII', 'CIV', 'OVI']
#mark=['v','o','x']
#colors=['fuchsia', 'blue', 'black']
#line=['-','.']


# Set the number of bins along each axis
numbins = 100
ewmin = -2.5

# Location of the master large ALL files
file_loc = '/home/matrix3/jrvander/sebass_gals/dwarfs/ALLfiles/masters/'
#ile_loc = '/home/jacob/matrix/sebass_gals/dwarfs/ALLfiles/masters/'


histos = []
xed = []
yed = []

# Loop over the ions
for i in range(0,len(galID_list)):

    if bulk!=0:
        # Create the four subplots
        fig,((p11,p12,p13), (p21,p22,p23), 
             (p31,p32,p33), (p41,p42,p43)) = plt.subplots(4,3,figsize=(10.2,10.8))
        plot_list = [p11, p21, p31, p41, p12, p22, p32, p42, p13, p23, p33, p43]

    # Loop over feedback models
    for ion in ion_list:    
        galID = galID_list[i]
    
        # Read in the master ALL file
        abs_file = galID+'.'+ion+'.ALL.sysabs.large.master'
        EW, lognH, logT, impact = np.loadtxt(file_loc+abs_file, skiprows=1, usecols=(5, 7, 8,  22), unpack=True)

        

        # Take the log of the EW array if do_logEW is high
        if do_logEW==1:
            EW = np.log10(EW)
            
            # Replace all NaNs in EW with zeros
            EW[np.isinf(abs(EW))] = -5.0


        # Bin the data
        H, xedges, yedges = np.histogram2d( impact, EW, bins=numbins, normed=normalize)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        

        # Convert the histogram to covering fractions
        coveringFrac = np.copy(H)

        # Loop through the impact parameter bins to determine 
        # how many lines of sight are in each bin
        numLOS = np.zeros(numbins)
        for j in range(0,numbins):
            minimp, maximp = xedges[j], xedges[j+1]
            count = 0.
            for D in impact:
                if D>=minimp and D<maximp:
                    count += 1
            numLOS[j] = count

        # Loop over impact bins
        for j in range(0,numbins):
            nLOS = numLOS[j]
            
            # Loop over EW bins
            for k in range(0,numbins):
                hits = 0

                # Loop over all bins with EW > the EW of H[:,k]
                for l in range(k, numbins):
                    hits += H[j,l]

                # Compute the covering fraction
                coveringFrac[j,k] = hits/nLOS

        H = coveringFrac

        # Rotate and filp the histogram
        H = np.rot90(H)
        H = np.flipud(H)
        
        # Mask the bins where the count is zero
        Hmasked = np.ma.masked_where(H==0,H)
        
        # Take the log of the count if do_logCount is high
        if do_logCount==1:
            Hmasked = np.log10(Hmasked)

         
        # Plot and save the figure
        if bulk==0:
            fig = plt.figure()
            plt.pcolormesh(xedges,yedges,Hmasked)
            plt.xlabel(' $Impact Parameter [$R_{vir}$] ')
            plt.ylabel(' EW [$\AA$] ')
            plt.xlim([0, 1.5])
            plt.ylim([ewmin,8])
            plt.title(feedback_list[i] +', '+ ion)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('$\log$ (Counts)')
            if normalize==0:
                fig.savefig(file_loc+'ewImp_histo.'+galID+'.'+ion+'.pdf')
                fig.savefig(file_loc+'ewImp_histo.'+galID+'.'+ion+'.eps')
                fig.savefig(file_loc+'ewImp_histo.'+galID+'.'+ion+'.jpg')
            else:
                fig.savefig(file_loc+'ewImp_histo.'+galID+'.'+ion+'.normalized.pdf')
                fig.savefig(file_loc+'ewImp_histo.'+galID+'.'+ion+'.normalized.eps')
                fig.savefig(file_loc+'ewImp_histo.'+galID+'.'+ion+'.normalized.jpg')

            plt.cla()
            plt.clf()

        else:
            xed.append(xedges)
            yed.append(yedges)
            histos.append(Hmasked)
            continue

            ax = plot_list[ion_list.index(ion)]
            mesh = ax.pcolormesh(xedges,yedges,Hmasked)
            ax.set_xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
            ax.set_ylabel(' $\log$ (T) [K] ')
            ax.set_xlim([0, 1.5])
            ax.set_ylim([ewmin,8])
            ax.text(-1, 7, ion)
            cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
            cbar.ax.set_label('$\log$ (Counts)')


# Plot if bulk is high
labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']

labels = ['(a)', '(d)', '(g)', '(j)', 
          '(b)', '(e)', '(h)', '(k)', 
          '(c)', '(f)', '(i)', '(l)']
ewmax = [1.5, 0.5, 0.0, 0.0,
         1.5, 0.5, 0.0, 0.0, 
         1.5, 0.5, 0.0, 0.0]


if bulk!=0:
    for i in range(0,len(plot_list)):
        ax = plot_list[i]
        H = histos[i]
#        H = np.log10(H)
        xedges = xed[i]
        yedges = yed[i]
        
        mesh = ax.pcolormesh(xedges,yedges,H)
        ax.set_xlim([0.0, 1.5])
        
        if do_logEW==1:
            ax.set_ylim([-2.5, ewmax[i]])

#        ax.text(-1, 7, labels[i])
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)


        if i==0:
            ax.set_title('dwSN', size=12)
            ax.set_ylabel('$\log$ EW$_{HI}$ [$\AA$] ', multialignment='center')
        elif i==1:
            ax.set_ylabel('$\log$ EW$_{MgII}$ [$\AA$] ', multialignment='center')
        elif i==2:
            ax.set_ylabel('$\log$ EW$_{CIV}$ [$\AA$] ', multialignment='center')
        elif i==3:
            ax.set_ylabel('$\log$ EW$_{OVI}$ [$\AA$] ', multialignment='center')
            ax.set_xlabel(' $D/R_{vir}$ ')
        elif i==4:
            ax.set_title('dwALL_1', size=12)
        elif i==7:
            ax.set_xlabel(' $D/R_{vir}$ ')
        elif i==8:
            ax.set_title('dwALL_8', size=12)
        elif i==11:
            ax.set_xlabel(' $D/R_{vir}$ ')

        if do_logCount==1:
            cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True, 
                                format=ticker.FuncFormatter(fmt))
        else:
            cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True, 
                                format=ticker.FuncFormatter(fmt2))
#        cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
#        cbar.ax.set_label('$\log$ (Counts)')
#        cbar.ax.tick_params(labelsize=10)

        cl = cbar.ax.get_yticklabels()

    dropx = [p12, p22, p32, p13, p23, p33, p11, p21, p31]
    plt.setp([a.get_xticklabels() for a in dropx],visible=False)

    dropy = [p12, p22, p32, p13, p23, p33, p42, p43]
    plt.setp([a.get_yticklabels() for a in dropy],visible=False)

    # Label the rows
#    fig.text(0.05, 0.85, 'HI', rotation='vertical')
#    fig.text(0.05, 0.63, 'MgII', rotation='vertical')
#    fig.text(0.05, 0.39, 'CIV', rotation='vertical')
#    fig.text(0.05, 0.15, 'OVI', rotation='vertical')

    plt.tight_layout()
    s = file_loc+'ewImp_cover.master_bulk.eps'
#    plt.savefig(s, bbox_inches='tight')

    if do_logEW==1 and do_logCount==1:
        s = file_loc+'ewImp_cover.master_bulk.logEW.logCount.pdf'
    elif do_logEW==0 and do_logCount==1:
        s = file_loc+'ewImp_cover.master_bulk.linEW.logCount.{0:d}.pdf'.format(numbins)
    elif do_logEW==1 and do_logCount==0: 
        s = file_loc+'ewImp_cover.master_bulk.logEW.linCount.{0:d}.pdf'.format(numbins)
    else: 
        s = file_loc+'ewImp_cover.master_bulk.linEW.linCount.pdf'



    plt.savefig(s, bbox_inches='tight')
