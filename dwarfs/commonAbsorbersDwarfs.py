from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import itertools as it
import scipy.stats as st
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator



kpc2km = 3.086e16 
pc2cm = 3.086e18 
kpc2cm = 3.086e21
mH = 1.6737e-24 
mSun = 1.989e33 
s2yr = 3.154e7 


def fmt(x,pos):
    return '{0:.1f}'.format(x)

def mkHist(ax,x,y,z,stat,xlabel,ylabel,cbarax,ion):
    numbins = 100
    #stat = 'count'
    if 'z' in ylabel:
        binrange = [[-3,3],[0,6]]
    else:
        binrange = [[-8,1],[2,8]]
    h,xedges,yedges,binnumber = st.binned_statistic_2d(x,y,z,
                              statistic=stat, bins=numbins,
                                 range=binrange)
    h = np.rot90(h)
    h = np.flipud(h)
    h[np.isnan(h)] = 0.0
    h = np.ma.masked_where(h==0,h)
    if stat=='sum':
        h = h/6.
    h = np.log10(h)
 
    if stat=='sum':
        if ion=='HI':
            vmin,vmax = 0,7
        else:
            vmin,vmax = 0,6
    elif stat=='mean':
        vmin,vmax = 0,5

    mesh = ax.pcolormesh(xedges,yedges,h,cmap='viridis',vmin=vmin,vmax=vmax)
    ax.set_xlim(binrange[0])
    ax.set_ylim(binrange[1])
    ax.set_xlabel(xlabel,fontsize='x-large')
    ax.set_ylabel(ylabel,fontsize='x-large')
    
    print(h.max(),h.min())
    #return h,xedges,yedges,binnumber
    return mesh

onMac = 0
if onMac==1:

    # ## Data Read In On Mac
    loc = '/Users/jacob/research/dwarfs/abscells/'

    filename = '{0:s}.*.{1:s}.abs_cells.dat'
    ions = 'HI MgII CIV OVI'.split()
    galIDs = 'D9o2 D9q D9m4a'.split()

    header = 'los impact cellID cellz logN b r nH t size SNII SNIa Zmet'.split()
    fullheader = 'los impact cellID cellz logN b r nH t size SNII SNIa Zmet galID ion boxz'.split()

    df = pd.DataFrame(columns=fullheader)
    for ion in ions:
        for galID in galIDs:
            print(ion,galID)
            fname = loc+filename.format(galID,ion)
            files = glob.glob(fname)
            for f in files:
                fparts = f.split('.')
                zfloat = float(fparts[1])+float(fparts[2])/1000.
                z = '{0:.3f}'.format(zfloat)
                d = np.loadtxt(f,skiprows=1)
                try:
                    d = pd.DataFrame(d,columns=header)
                    d['galID'] = galID
                    d['ion'] = ion
                    d['boxz'] = z
                    df = df.append(d)
                except ValueError:
                    continue
            
else: 
    # ## Data Read In On Linux 
    loc = '/home/jacob/research/dwarfs/abscells/'
    loc = '/mnt/cluster/abs/cgm/dwarfs/'

    filename = '{0:s}.{1:s}.bulk_abscells.dat'
    ions = 'HI MgII CIV OVI'.split()
    galIDs = 'D9o2 D9q D9m4a'.split()

    header = 'los impact cellID cellz logN b r nH t size SNII SNIa Zmet boxz'.split()
    fullheader = 'los impact cellID cellz logN b r nH t size SNII SNIa Zmet galID ion boxz'.split()

    df = pd.DataFrame(columns=fullheader)
    for ion in ions:
        for galID in galIDs:
            print(ion,galID)
            fname = loc+filename.format(galID,ion)
            d = np.loadtxt(fname,skiprows=1)
            try:
                d = pd.DataFrame(d,columns=header)
                d['galID'] = galID
                d['ion'] = ion
                df = df.append(d)
            except ValueError:
                continue


# In[ ]:

groups = df.groupby(['galID','boxz'])

ionsCombo = it.combinations(ions,2)
common1 = pd.DataFrame(columns=fullheader)
common2 = pd.DataFrame(columns=fullheader)
common3 = pd.DataFrame(columns=fullheader)
totals = np.zeros((len(ions),len(ions)))
totals = pd.DataFrame(totals,columns=ions,index=ions)
for ion1,ion2 in ionsCombo:
    print(ion1,ion2)
    for (gal,a),group in groups:
        cells1 = group['cellID'][group['ion']==ion1]
        cells2 = group['cellID'][group['ion']==ion2]
        common = set(cells1) & set(cells2)
        commondf = group[(group['cellID'].isin(common)) & (group['ion']==ion1)]
        totals[ion1].ix[ion2] += len(commondf)
        totals[ion2].ix[ion1] += len(commondf)
        totals[ion1].ix[ion1] += len(cells1)/3.
        commondf['ion'] = '{0:s} {1:s}'.format(ion1,ion2)
        if gal=='D9o2':
            common1 = common1.append(commondf)
        elif gal=='D9q':
            common2 = common2.append(commondf)
        else:
            common3 = common3.append(commondf)


# In[ ]:

titles = 'D9SN D9ALL\_1 D9ALL\_8'.split()
ionsCombo = it.combinations(ions,2)
for ion1,ion2 in ionsCombo:
    print(ion1,ion2)
    fig,axes = plt.subplots(1,3,figsize=(15,5),sharey='row',squeeze=True)
    cbarax = fig.add_axes([.94,.2,.03,.6])
    s = '{0:s} {1:s}'.format(ion1,ion2)
    for ax,com,title in zip(axes,[common1,common2,common3],titles):
        df = com[com['ion']==s]
        df['mass'] = 10**df['nH']*mH*(df['size']*kpc2cm)**3 / mSun
        mesh = mkHist(ax,df['nH'],df['t'],df['mass'],'sum','','',cbarax)
        ax.annotate(title,(-2,7),xycoords='data',fontsize='x-large')
        ax.set_xlabel('log( $n_H$ [cm$^{-3}$])',fontsize='x-large')
    axes[0].set_ylabel('log( T [K] )',fontsize='x-large')
    axes[0].xaxis.set_major_locator(MaxNLocator(prune='upper'))          
    axes[1].xaxis.set_major_locator(MaxNLocator(prune='upper'))          
    cbar = plt.colorbar(mesh,cax=cbarax,use_gridspec=True,
                        format=ticker.FuncFormatter(fmt))
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('log (M$_g$ [M$_{\odot}$])',rotation=270,fontsize='x-large')
    fig.subplots_adjust(wspace=0)
    s = 'commonAbs_{0:s}_{1:s}.pdf'.format(ion1,ion2)
    fig.savefig(s,bbox_inches='tight',dpi=300)


# ## Full Phase 

# In[ ]:

runs = 'D9SN D9ALL\_1 D9ALL\_8'.split()
fig,axes = plt.subplots(4,3,figsize=(15,20),sharey='row',sharex='col')
cbarax1 = fig.add_axes([.94,.72,.03,.17])
cbarax2 = fig.add_axes([.94,.53,.03,.17])
cbarax3 = fig.add_axes([.94,.33,.03,.17])
cbarax4 = fig.add_axes([.94,.13,.03,.17])
cbaraxes = [cbarax1,cbarax2,cbarax3,cbarax4]
meshes = []

for j,ion in enumerate(ions):
    for i,galID in enumerate(galIDs):
        index = (df['ion']==ion) & (df['galID']==galID)
        d = df[index]
        d['mass'] = 10**d['nH']*mH*(d['size']*kpc2cm)**3 / mSun
        ax = axes[j,i]
        mesh = mkHist(ax,d['nH'],d['t'],d['mass'],'sum','','',cbarax1,ion)
        s = '{0:s},{1:s}'.format(galID,ion)
        #ax.annotate(s,(-2,7),xycoords='data')
        meshes.append(mesh)
for i,cbarax in enumerate(cbaraxes):
    mesh = meshes[3*(i+1)-1]
    cbar = plt.colorbar(mesh,cax=cbarax,use_gridspec=True,
                    format=ticker.FuncFormatter(fmt))
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('log (M$_g$ [M$_{\odot}$])',rotation=270,fontsize='x-large')    

    #cbar.ax.set_ylabel('{0:d}'.format(i),rotation=270,fontsize='x-large')    


for ax in axes[3,:]:
    ax.set_xlabel('log( $n_H$ [cm$^{-3}$])',fontsize='x-large')
for ax,ion in zip(axes[:,0],ions):
    ax.set_ylabel('{0:s}\nlog( T [K] )'.format(ion),fontsize='x-large') 
for ax,run in zip(axes[0,:],runs):
    ax.set_title(run,fontsize='x-large')
    
axes[3,0].xaxis.set_major_locator(MaxNLocator(prune='upper'))
axes[3,1].xaxis.set_major_locator(MaxNLocator(prune='upper'))

#axes[3,0].yaxis.set_major_locator(MaxNLocator(prune='upper'))
#axes[2,0].yaxis.set_major_locator(MaxNLocator(prune='upper'))
#axes[1,0].yaxis.set_major_locator(MaxNLocator(prune='upper'))

axes[1,0].set_yticks([2,3,4,5,6,7])
axes[2,0].set_yticks([2,3,4,5,6,7])
axes[3,0].set_yticks([2,3,4,5,6,7])
    
fig.subplots_adjust(wspace=0,hspace=0)
s = 'dwarfFullPhase.pdf'
fig.savefig(s,bbox_inches='tight',dpi=100)
plt.close(fig)


# ## Probabilities 

probs = np.zeros((len(ions),len(ions)))
probs = pd.DataFrame(probs,columns=ions,index=ions)
for givenIon in ions:
    for testIon in ions:
        if testIon!=givenIon:
            numCommon = 0.
            numTotal = 0.
                
            for galID in galIDs:
                d = df[df['galID']==galID]
                expns = d['boxz'].unique()                

                for expn in expns:
                    box = d[d['boxz']==expn]
                    givenIndex = (box['ion']==givenIon)
                    numTotal += givenIndex.sum()
                    givenCells = box['cellID'][givenIndex]
                    testIndex = (box['ion']==testIon)
                    testCells = box['cellID'][testIndex]
                    numCommon += len(set(givenCells) & set(testCells))
            frac = numCommon/numTotal
            print('Given {0:s}, fraction with {1:s} = {2:.3f}'.format(givenIon,testIon,frac))
            probs[testIon].ix[givenIon] = frac


# In[ ]:

probs.to_csv('commonAbsDwarfsProps.dat',sep='\t')


# ## Summary File

# In[ ]:

outname = 'dwarfCommonAbsSummary.h5'
header = 'galID expn cellID impact r nH t HI MgII CIV OVI'.split()
common = pd.DataFrame(columns=header)
for galID in galIDs:
    d = df[df['galID']==galID]
    expns = d['boxz'].unique()
    print(galID)
    for expn in expns:
        print(expn,type(expn))
        box = d[d['boxz']==expn]
        cells = box['cellID'].unique()
        print(len(box))
        for cell in cells:
            absData = pd.Series(index=header)
            index = box['cellID']==cell
            absorbing = box[box['cellID']==cell]
            absIons = absorbing['ion']
            for ion in absIons:
                absData[ion] = 1
            absData['galID'] = galID
            absData['expn'] = expn
            absData['cellID'] = cell
            absData['impact'] = absorbing['impact']
            absData['r'] = absorbing['r']
            absData['nH'] = absorbing['nH']
            absData['t'] = absorbing['t']
            absData.fillna(0,inplace=True)
            common = common.append(absData,ignore_index=True)

# In[ ]:

common.to_hdf('dwarfCommonAbsSummary.h5','data',mode='w')


