
# coding: utf-8

# In[1]:

from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().magic(u'matplotlib inline')
pd.options.mode.chained_assignment = None
sns.set(style='white')


# ### Individual EW vs Impact

# In[168]:

loc = '/home/jacob/research/code/analysis/bulkVelas/'
filename = 'vela2b-{0:d}_ewImpact.h5'
galNums = range(21,30)
ions = 'HI MgII CIV OVI'.split()


# In[169]:

for i,galNum in enumerate(galNums):
    fig,axes =  plt.subplots(1,4,figsize=(20,5),sharey='row')
    for ion,ax in zip(ions,axes):
        fname = loc+filename.format(galNum)
        df = pd.read_hdf(fname,'data')
        hi = df.xs(ion,level=1,axis=1).transpose()
        hi = np.log10(hi.mask(hi==0))
        sns.heatmap(hi,ax=ax,vmin=-2,vmax=1,cmap='viridis',annot=False,fmt='0.1f',cbar=(ion=='OVI'))
        ax.set_title(ion)
        for spine in ['left','right','top','bottom']:
            ax.spines[spine].set_visible(True)
            ax.spines[spine].set_linewidth(0.75)
            ax.spines[spine].set_color('#262626')

    fig.subplots_adjust(wspace=0)
    s = '{0:s}/ewDist_vela2b-{1:d}.png'.format(loc,galNum)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


# ### Mean EW vs Impact

# In[171]:

fig,axes = plt.subplots(1,4,figsize=(20,5),sharey='row')
cbarax = fig.add_axes([.91,.125,.01,.775])
for ion,ax in zip(ions,axes):
    for i,galNum in enumerate(galNums):
        fname = loc+filename.format(galNum)
        df = pd.read_hdf(fname,'data')
        hi = df.xs(ion,level=1,axis=1).transpose()
        
        expns = [float(a.split('a')[-1])/1000. for a in hi.index]
        hi.index = ['{0:.2f}'.format(1./a-1) for a in expns]
        
        if i==0:
            totals = hi
        else:
            totals = totals.add(hi,fill_value=0)
    totals /= float(i+1)
    m1 = np.log10(totals.mask(totals==0).fillna(100).values.min())
    m2 = np.log10(totals.fillna(-4).values.max())
    print(m1,m2)

    totals = np.log10(totals.mask(totals==0))
    sns.heatmap(totals[2:],ax=ax,vmin=-1,vmax=1,
                cmap='viridis',
                annot=False,fmt='0.1f',
                cbar_ax=cbarax,
                yticklabels=[a if i%2==0 else '' for i,a in enumerate(totals.index)])
                #cbar=(ion=='OVI'))
    for spine in ['left','right','top','bottom']:
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_linewidth(0.75)
        ax.spines[spine].set_color('#262626')
    ax.invert_yaxis()
    ax.set_title(ion)
    ax.set_xlabel('Impact Parameter [Rvir]')
axes[0].set_ylabel('Redshift')
#axes[0].set_yticks(rotation=0)
plt.setp( axes[0].yaxis.get_majorticklabels(), rotation=0 )
cbarax.get_yaxis().labelpad = 20
cbarax.set_ylabel(r'log (EW [\AA])',rotation=270)#,fontsize='large')    

#axes[0].set_yticks([a if i%2==0 else '' for i,a in enumerate(totals.index)])
#axes[0].set_yticks([0]*len(hi))
#axes[0].set_yticks(zticks)
print(zticks)
print(len(zticks),len(totals))
fig.subplots_adjust(wspace=0.05)
s = '{0:s}/ewDist_vela2b-Mean.png'.format(loc)
fig.savefig(s,bbox_inches='tight',dpi=300)


# # Covering Fraction

# In[141]:

loc = '/home/jacob/research/code/analysis/bulkVelas/'
filename = 'vela2b-{0:d}_covering.h5'
galNums = range(21,30)
ions = 'HI MgII CIV OVI'.split()


# In[148]:

for i,galNum in enumerate(galNums):
    fig,axes =  plt.subplots(1,4,figsize=(20,5),sharey='row')
    cbarax = fig.add_axes([.91,.125,.03,.775])
    for ion,ax in zip(ions,axes):
        fname = loc+filename.format(galNum)
        df = pd.read_hdf(fname,'data')
        hi = df.xs(ion,level=1,axis=1).transpose()
        hi = hi.mask(hi==0)
        sns.heatmap(hi,ax=ax,vmin=0,vmax=1,cmap='viridis',
                    annot=False,fmt='0.1f',
                    cbar_ax=cbarax,
                    yticklabels=[a if i%2==0 else '' for i,a in enumerate(totals.index)])
        #ax.invert_yaxis()
        ax.set_title(ion)
        ax.set_xlabel('Impact Parameter [Rvir]')
    axes[0].set_ylabel('Redshift')
    plt.setp( axes[0].yaxis.get_majorticklabels(), rotation=0 )
    cbarax.get_yaxis().labelpad = 20
    cbarax.set_ylabel('$C_f$',rotation=270)#,fontsize='large')  
    fig.subplots_adjust(wspace=0.5)
    s = '{0:s}/covering_vela2b-{1:d}.png'.format(loc,galNum)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


# In[162]:

fig,axes = plt.subplots(1,4,figsize=(20,5),sharey='row')
cbarax = fig.add_axes([.91,.125,.01,.775])
for ion,ax in zip(ions,axes):
    for i,galNum in enumerate(galNums):
        fname = loc+filename.format(galNum)
        df = pd.read_hdf(fname,'data')
        hi = df.xs(ion,level=1,axis=1).transpose()
        
        expns = [float(a.split('a')[-1])/1000. for a in hi.index]
        hi.index = ['{0:.2f}'.format(1./a-1) for a in expns]
        if i==0:
            totals = hi
        else:
            totals = totals.add(hi,fill_value=0)
    totals /= float(i+1)

    totals = totals.mask(totals==0)
    print(totals.index)
    sns.heatmap(totals[:-4],ax=ax,vmin=0,vmax=1,
                cmap='viridis',
                annot=False,fmt='0.1f',
                cbar_ax=cbarax,
                yticklabels=[a if i%2==0 else '' for i,a in enumerate(totals.index[:-4])])
                #cbar=(ion=='OVI'))
    ax.set_title(ion)
    #ax.invert_yaxis()
    ax.set_xlabel('Impact Parameter [Rvir]')
    
    for spine in ['left','right','top','bottom']:
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_linewidth(0.75)
        ax.spines[spine].set_color('#262626')
axes[0].set_ylabel('Redshift')
plt.setp( axes[0].yaxis.get_majorticklabels(), rotation=0 )
cbarax.get_yaxis().labelpad = 20
cbarax.set_ylabel('$C_f$',rotation=270)#,fontsize='large')    
fig.subplots_adjust(wspace=0.05)
s = '{0:s}/covering_vela2b-Mean.png'.format(loc)
fig.savefig(s,bbox_inches='tight',dpi=300)


# In[ ]:




# In[51]:

t.fillna(-4).values.max()


# In[77]:

expns = [float(a.split('a')[-1])/1000. for a in totals.index]
totals.index = ['{0:.2f}'.format(1./a-1) for a in expns]


# In[78]:

totals


# In[102]:

aticks = np.arange(0.200,0.500,0.02)
zvals = ['{0:.2f}'.format(1./a-1) for a in aticks]


# In[103]:

zblanks = ['']*len(zticks)


# In[104]:

zticks = [val for pair in zip(zvals, zblanks) for val in pair]


# In[106]:

zticks


# In[ ]:



