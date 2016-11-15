
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools as it


def plotlines(ax,name):
    if 'vr' in name:
        lower,upper = 0,1
    else:
        lower,upper = 0,5
    ax.vlines(0.29,lower,upper,colors='r',linestyles='dashed',linewidth=1)
    ax.vlines(0.32,lower,upper,colors='r',linestyles='dashed',linewidth=1)
    ax.vlines(0.35,lower,upper,colors='g',linestyles='dashed',linewidth=1)
    ax.vlines(0.45,lower,upper,colors='g',linestyles='dashed',linewidth=1)

tempLabels = 'cool warm hot'.split()

# In[3]:

dataloc = '/home/jacob/research/code/analysis/inflows/'  
filebase = 'vradFrac r'.split()
stats = 'mean std'.split()
stores,names = [],[]
for combos in it.product(filebase,stats):
    filename = '{0:s}_{1:s}.h5'.format(combos[0],combos[1])
    store = pd.HDFStore(filename)
    stores.append(store)
    names.append('{0:s}{1:s}'.format(combos[0],combos[1].title()))
    

#vrMean = pd.read_hdf(dataloc+'vradFraction_mean.h5','data')
#vrStd = pd.read_hdf(dataloc+'vradFraction_std.h5','data')
#rMean = pd.read_hdf(dataloc+'r_mean.h5','data')
#rStd = pd.read_hdf(dataloc+'r_std.h5','data')
#frames = [vrMean,vrStd,rMean,rStd]
#names = 'vrMean vrStd rMean rStd'.split()


# In[51]:

for i,tempLabel in enumerate(tempLabels):
    fig,axes = plt.subplots(2,2,figsize=(10,10))
    axes = axes.flatten()
    m = 's ^ o x * D'.split()
    for j,ax in enumerate(axes):
        df = stores[j][tempLabel]
        for k,ind in enumerate(df.columns):
            line = ax.plot(np.abs(df[ind]), linestyle='solid',
                    marker=m[k],label='${0:s}$'.format(ind))
        plotlines(ax,names[j])
        ax.set_title(names[j])
    handles,labels = ax.get_legend_handles_labels()
    #plt.legend(handles, labels,bbox_to_anchor=(1,2.250),ncol=6,
    #          borderaxespad=0.,fontsize='large')

    ax.legend(handles,labels,loc='best')
    #fig.legend(lines,labels,'right')#bbox_to_anchor=(1,0.5),ncol=1)
    fig.tight_layout()
    fig.suptitle(tempLabel,fontsize='x-large')
        
    s = 'velFrac_{0:s}.png'.format(tempLabel)
    fig.savefig(s,bbox_inches='tight',dpi=300)

for store in stores:
    store.close()


