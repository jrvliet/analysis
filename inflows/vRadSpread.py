
# coding: utf-8

# In[1]:

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None
import seaborn as sns


# In[9]:

codeloc = '/home/jacob/research/code/analysis/inflows/'
dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
expns = np.arange(0.200,0.550,0.01)


# ### Get virial radius

# In[10]:

virialFile = codeloc+'vela2b-27_virialProperties.h5'
virial = pd.read_hdf(virialFile,'data')


# ### Define bins

# In[11]:

lowestT,highestT = 3.5,6.5
nTempBins = 3
tempBins = np.linspace(lowestT,highestT,nTempBins+1)
loN,hiN = -5.5,-2.5
numRBins = 100
numRBins = 100
rLabels = range(numRBins)
rBins = np.linspace(0,6,numRBins)

# Dataframes to hold final results
header = 'a vrMean_Mean vrMean_Median vrStd_Mean vrStd_Median'.split()
index = range(len(expns))
coolSummary = pd.DataFrame(columns=header,index=index)
warmSummary = pd.DataFrame(columns=header,index=index)
hotSummary = pd.DataFrame(columns=header,index=index)
summaries = [coolSummary, warmSummary, hotSummary]


# In[12]:

store = pd.HDFStore(codeloc+'vRadSpreadEvolution.h5',mode='w')
for i,a in enumerate(expns):
    
    print(a)
    # Get rvir
    rvir = virial[np.isclose(virial['a'],a)]['rvir'].values[0]
    
    # Read in data
    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname,'data')

    denseBins = (df['density']>10**loN) & (df['density']<10**hiN)
    spaceBins = (df['theta']<80)

    tempHeader = 'cool warm hot'.split()
    fieldHeader = 'rModMean vrMean vrStd vrIQR numCells'.split()
    headers = [tempHeader,fieldHeader]
    header = pd.MultiIndex.from_product(headers)#,names=['temperature','quantity'])
    data = np.zeros((numRBins,len(header)))
    vRadSpread = pd.DataFrame(data,columns=header,index=rLabels)

    for j in range(nTempBins):
        loT = tempBins[j]
        hiT = tempBins[j+1]
        tempLabel = tempHeader[j]
        print('\t{0:s}'.format(tempLabel))
        temperatureBins = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
        fil = df[denseBins & spaceBins & temperatureBins]

        # Get r and vr
        fil['r'] = np.sqrt(fil['xRot']**2 + fil['yRot']**2 + fil['zRot']**2)
        fil['vr'] = (fil['xRot']*fil['vxRot'] + fil['yRot']*fil['vyRot'] + 
                    fil['zRot']*fil['vzRot']) / fil['r']
        fil['rMod'] = fil['r']/rvir

        # Bin data by rMod
        fil['rBin'] = pd.cut(fil['rMod'],numRBins,labels=rLabels)

        vrMean = np.zeros(numRBins)
        rModMean = np.zeros_like(vrMean)
        for rLabel in rLabels:
            vrMean[rLabel] = fil[fil['rBin']==rLabel]['vr'].mean()
            rModMean[rLabel] = fil[fil['rBin']==rLabel]['rMod'].mean()
            vRadSpread[tempLabel,'vrMean'].iloc[rLabel] = fil[fil['rBin']==rLabel]['vr'].mean()
            vRadSpread[tempLabel,'rModMean'].iloc[rLabel] = fil[fil['rBin']==rLabel]['rMod'].mean()

        # Add fields describing the mean r and vr for the bin that the cell 
        # belongs to, then subtract off the mean
        fil['rModMean'] = rModMean[fil['rBin']]
        fil['vrMean'] = vrMean[fil['rBin']]
        fil['vrDeTrend'] = fil['vr']-fil['vrMean']

        sigma = np.zeros_like(vrMean)
        iqr = np.zeros_like(vrMean)
        numCells = np.zeros_like(vrMean)
        for rLabel in rLabels:
            cells = (fil['rBin']==rLabel)
            iqrs = fil['vrDeTrend'][cells].quantile([0.25,0.75])
            vRadSpread[tempLabel,'vrStd'].iloc[rLabel] = fil['vrDeTrend'][cells].std()
            vRadSpread[tempLabel,'vrIQR'].iloc[rLabel] = iqrs[0.75]-iqrs[0.25]
            vRadSpread[tempLabel,'numCells'].iloc[rLabel] = cells.sum()

        # Store the result
        # Determine the summary properties
        summary = summaries[j]
        summary['a'].iloc[i] = a
        summary['vrMean_Mean'].iloc[i] = vRadSpread[tempLabel,'vrMean'].mean()
        summary['vrMean_Median'].iloc[i] = vRadSpread[tempLabel,'vrMean'].median()
        summary['vrStd_Mean'].iloc[i] = vRadSpread[tempLabel,'vrStd'].mean()
        summary['vrStd_Median'].iloc[i] = vRadSpread[tempLabel,'vrStd'].median()
        
    # Store the result
    storeLabel = 'a{0:d}'.format(int(a*1000))
    store[storeLabel] = vRadSpread
store.close()

store = pd.HDFStore(codeloc+'vRadSummaries.h5',mode='w')
for i,tempLabel in enumerate(tempHeader):
    store[tempLabel] = summaries[i]
store.close()


# In[ ]:

store = pd.HDFStore(codeloc+'vRadSummaries.h5')
fig,(ax1,ax2) = plt.subplots()
for s in store:
    df = store[s]
    ax1.plot(df['a'],df['vrMean_Mean'],label='{0:s}, Mean'.format(s[1:]))
    ax1.plot(df['a'],df['vrMean_Median'],label='{0:s}, Median'.format(s[1:]))
    ax2.plot(df['a'],df['vrStd_Mean'],label='{0:s}, Mean'.format(s[1:]))
    ax2.plot(df['a'],df['vrStd_Median'],label='{0:s}, Median'.format(s[1:]))
store.close()
ax1.set_xlabel('Expn')
ax2.set_xlabel('Expn')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax1.set_ylabel('Vr Mean')
ax2.set_ylabel('Vr Std')
fig.tight_layout()
s = 'vRadSummariesEvolution.png'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)


