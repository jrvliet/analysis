
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
import seaborn as sns


# In[18]:

loc = '/Users/jacob/research/code/analysis/inflows/'
filename = 'vRadSpreadEvolution.h5'
fname = loc+filename
store = pd.HDFStore(fname)


# In[25]:

store['a200'].head()


# In[49]:

tempLabels = 'cool warm hot'.split()
colors = 'blue green red'.split()
mins = [-250,-175,-125]
maxs = [125,275,375]

for s in store:
    df = store[s]
    a = float(s.split('a')[-1])/1000.
    print(a)
    fig,axes = plt.subplots(1,3,figsize=(15,5))
    for ax,tempLabel,color,vMin,vMax in zip(axes,tempLabels,colors,mins,maxs):
        label = '{0:s} Mean'.format(tempLabel.capitalize())
        x = df[tempLabel,'rModMean']
        y = df[tempLabel,'vrMean']
        ax.plot(x,y,color=color,label=label)
        ax.legend(loc='best')
        ax.set_xlabel('r [Rvir]')
        ax.set_ylabel('Vr [km/s]')
        ax.set_ylim([-250,400])
        ax.set_xlim([0,6])
        ax.hlines(0,0,6,colors='black',linewidth=2)
    axes[1].set_title('a={0:.3f}, z={1:.3f}'.format(a,1./a-1))
    name = 'vrMeanProfile_uniformAxis_a{0:.3f}.png'.format(a)
    fig.tight_layout()
    fig.savefig(name,bbox_inches='tight',dpi=300)
    plt.close(fig)


# In[33]:

a = np.zeros(len(store))
vrCoolMean = np.zeros_like(a)
vrWarmMean = np.zeros_like(a)
vrHotMean = np.zeros_like(a)
vrCoolMed = np.zeros_like(a)
vrWarmMed = np.zeros_like(a)
vrHotMed = np.zeros_like(a)
vrStdCool = np.zeros_like(a)
vrStdWarm = np.zeros_like(a)
vrStdHot = np.zeros_like(a)
vrRatioCool = np.zeros_like(a)
vrRatioWarm = np.zeros_like(a)
vrRatioHot = np.zeros_like(a)
for i,s in enumerate(store):
    df = store[s]
    a[i] = float(s.split('a')[-1])/1000.
    vrCoolMean[i] = df['cool','vrMean'].mean()
    vrWarmMean[i] = df['warm','vrMean'].mean()
    vrHotMean[i] = df['hot','vrMean'].mean()
    
    vrCoolMed[i] = df['cool','vrMean'].median()
    vrWarmMed[i] = df['warm','vrMean'].median()
    vrHotMed[i] = df['hot','vrMean'].median()

    vrStdCool[i] = df['cool','vrStd'].mean()
    vrStdWarm[i] = df['warm','vrStd'].mean()
    vrStdHot[i] = df['hot','vrStd'].mean()
    
    vrRatioCool[i] = (df['cool','vrStd']/np.abs(df['cool','vrMean'])).mean()
    vrRatioWarm[i] = (df['warm','vrStd']/np.abs(df['warm','vrMean'])).mean()
    vrRatioHot[i] = (df['hot','vrStd']/np.abs(df['hot','vrMean'])).mean()


# In[35]:

fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
ax1.plot(a,vrCoolMean,label='Cool Mean',color='blue',linestyle='solid')
ax1.plot(a,vrWarmMean,label='Warm Mean',color='green',linestyle='solid')
ax1.plot(a,vrHotMean,label='Hot Mean',color='red',linestyle='solid')
ax1.plot(a,vrCoolMed,label='Cool Med',color='blue',linestyle='dashed')
ax1.plot(a,vrWarmMed,label='Warm Med',color='green',linestyle='dashed')
ax1.plot(a,vrHotMed,label='Hot Med',color='red',linestyle='dashed')
ax1.legend(loc='best')

ax2.plot(a,vrStdCool,label='Cool',color='blue')
ax2.plot(a,vrStdWarm,label='Warm',color='green')
ax2.plot(a,vrStdHot,label='Hot',color='red')
ax2.legend(loc='best')

ax3.plot(a,vrRatioCool,label='Cool',color='blue')
#ax3.plot(a,vrRatioWarm,label='Warm',color='green')
ax3.plot(a,vrRatioHot,label='Hot',color='red')
ax3.legend(loc='best')
fig.tight_layout()


# In[42]:

get_ipython().magic('pwd')


# In[ ]:



