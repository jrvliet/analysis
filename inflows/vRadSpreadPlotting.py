
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#get_ipython().magic('matplotlib inline')
#import seaborn as sns
import sys

#sns.set_style('white')

# In[18]:

loc = '/Users/jacob/research/code/analysis/inflows/'
loc = '/home/jacob/research/code/analysis/inflows/'
filename = 'vRadSpreadEvolution.h5'
fname = loc+filename
store = pd.HDFStore(fname)

tempLabels = 'cool warm hot'.split()
colors = 'blue green red'.split()
lines = 'solid dashed dotted'.split()
widths = '1 1 2'.split()
mins = [-250,-175,-125]
maxs = [125,275,375]

#for s in store:
#    df = store[s]
#    a = float(s.split('a')[-1])/1000.
#    print(a)
#    fig,ax = plt.subplots(1,1,figsize=(5,5))
#    for tempLabel,color,line,width in zip(tempLabels,colors,lines,widths):
#        label = '{0:s} Mean'.format(tempLabel.capitalize())
#        x = df[tempLabel,'rModMean']
#        y = df[tempLabel,'vrMean']
#        ax.plot(x,y,color=color,label=label,linestyle=line,linewidth=width)
#
#    # Apply labels
#    if a>=0.29 and a<=0.32:
#        s = 'Merging'
#        ax.text(4.5,250,s)
#    elif a>0.32 and a<=0.35:
#        s = 'Starburst'
#        ax.text(4.5,250,s)
#    ax.legend(loc='upper right',fontsize='small')
#    ax.set_xlabel('r [Rvir]')
#    ax.set_ylabel('Vr [km/s]')
#    ax.set_ylim([-250,400])
#    ax.set_xlim([0,6])
#    ax.hlines(0,0,6,colors='black',linewidth=2)
#    ax.set_title('a={0:.3f}, z={1:.3f}'.format(a,1./a-1))
#    name = 'vrMeanProfile_uniformAxis_a{0:.3f}.png'.format(a)
#    fig.tight_layout()
#    fig.savefig(name,bbox_inches='tight',dpi=300)
#    plt.close(fig)
#
# In[33]:


# Plot only starburst curves
print('\n\nStarburst Only...')
sbStart = '/a320'
sbEnd = '/a350'
colors = 'blue green red'.split()
styles = '- -- -. :'.split()

width = 8
height = 10
fig1,axes1 = plt.subplots(4,1,figsize=(width,height))
fig2,axes2 = plt.subplots(4,1,figsize=(width,height))
fig3,axes3 = plt.subplots(4,1,figsize=(width,height))
figs = [fig1,fig2,fig3]
axes = [axes1,axes2,axes3]

i = 0

coolVr = store[sbStart]['cool','vrMean'].to_frame()
warmVr = store[sbStart]['warm','vrMean'].to_frame()
hotVr = store[sbStart]['hot','vrMean'].to_frame()

coolr = store[sbStart]['cool','rModMean'].to_frame()
warmr = store[sbStart]['warm','rModMean'].to_frame()
hotr = store[sbStart]['hot','rModMean'].to_frame()

for s in store:
    if s>sbStart and s<=sbEnd:

        vr = store[s]['cool','vrMean'].to_frame()
        r = store[s]['cool','rModMean'].to_frame()
        coolVr = pd.concat([coolVr,vr],axis=1)
        coolr = pd.concat([coolr,r],axis=1)

        vr = store[s]['warm','vrMean'].to_frame()
        r = store[s]['warm','rModMean'].to_frame()
        warmVr = pd.concat([warmVr,vr],axis=1)
        warmr = pd.concat([warmr,r],axis=1)

        vr = store[s]['hot','vrMean'].to_frame()
        r = store[s]['hot','rModMean'].to_frame()
        hotVr = pd.concat([hotVr,vr],axis=1)
        hotr = pd.concat([hotr,r],axis=1)

        a = float(s.split('a')[-1])/1000.
        print(a)
        redshift = 'z={0:.2f}'.format(1./a - 1.)

fig,ax = plt.subplots(1,1,figsize=(5,5))
# Cool gas
x = coolr.mean(axis=1)
y = coolVr.mean(axis=1)
err = coolVr.std(axis=1)
print('Cool: {0:.6f}'.format(err.mean()))
lower = y-err
upper = y+err
ax.plot(x,y,color='blue',label='cool')
ax.fill_between(x,lower,upper,color='blue',alpha=0.25)

# Warm gas
x = warmr.mean(axis=1)
y = warmVr.mean(axis=1)
err = warmVr.std(axis=1)
print('Warm: {0:.6f}'.format(err.mean()))
lower = y-err
upper = y+err
ax.plot(x,y,color='green',label='warm')
ax.fill_between(x,lower,upper,color='green',alpha=0.25)

# Hot gas
x = hotr.mean(axis=1)
y = hotVr.mean(axis=1)
err = hotVr.std(axis=1)
print('Hot: {0:.6f}'.format(err.mean()))
lower = y-err
upper = y+err
ax.plot(x,y,color='red',label='hot')
ax.fill_between(x,lower,upper,color='red',alpha=0.25)

ax.set_ylabel('Vr [km/s]')
ax.set_ylim([-250,400])
ax.set_xlim([0,5])
ax.hlines(0,0,6,colors='black',linewidth=2)
ax.set_xlabel('r [Rvir]')
ax.legend(ncol=3,loc='upper center')
s = 'starburst_vrMeanProfile_mean.pdf'
fig.savefig(s,bbox_inches='tight',dpi=300)
plt.close(fig)



store.close()
sys.exit()















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



