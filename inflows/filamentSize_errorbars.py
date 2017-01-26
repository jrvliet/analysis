
# coding: utf-8

# In[2]:

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[3]:

loc = '/home/jacob/research/code/analysis/inflows/'
errFilename = 'filamentSize_errorbars.h5'
errFname = loc+errFilename
valFilename = 'massContainedSimple.h5'
valFname = loc+valFilename


# In[4]:

errors = pd.read_hdf(errFname,'data')
valStore = pd.HDFStore(valFname)


# In[5]:

valStore


# In[6]:

valStore['nLim=(-5.5,-4.5)'].head()


# In[86]:

errors = errors.apply(pd.to_numeric)


# In[8]:

df = pd.read_hdf(loc+'filamentMassSize_simple.h5','data')


# In[11]:

expns = []
for a in df.index:
    expns.append(float(a.split('a')[-1])/1000.)


# In[79]:

errors.head()


# In[13]:

fig1,ax1 = plt.subplots(1,1,figsize=(5,5))
fig2,ax2 = plt.subplots(1,1,figsize=(5,5))

z = [1./a-1 for a in expns]

datLabs = 'rare mid dense'.split()
errLabs = 'low mid high'.split()
labels = 'diffuse mid dense'.split()
colors = 'blue green red'.split()

for datLab,errLab,label,c in zip(datLabs,errLabs,labels,colors):
    ax1.errorbar(z,df['cool',datLab,'rho90kpc'],yerr=errors[errLab,'stdpkpc'],label=label,color=c)
    ax2.errorbar(z,df['cool',datLab,'rho90Rvir'],yerr=errors[errLab,'stdRvir'],label=label,color=c)

ax1.set_ylabel(r'$r_{90}$ [pkpc]')
ax2.set_ylabel(r'$r_{90}$ [Rvir]')

for ax in (ax1,ax2):
    ax.legend(loc='best')
    ax.set_xlabel('Redshift')
    ax.invert_xaxis()
    ax.set_xlim([4,1])
    
fig1.savefig(loc+'filamentSizeErrors_pkpc.pdf',bbox_inches='tight',dpi=300)
fig2.savefig(loc+'filamentSizeErrors_rvir.pdf',bbox_inches='tight',dpi=300)
plt.close(fig1)
plt.close(fig2)


# In[92]:

fig1,ax1 = plt.subplots(1,1,figsize=(5,5))
fig2,ax2 = plt.subplots(1,1,figsize=(5,5))

z = [1./a-1 for a in expns]

datLabs = 'rare mid dense'.split()
errLabs = 'low mid high'.split()
labels = 'diffuse mid dense'.split()
colors = 'blue green red'.split()
print(np.isfinite(z).sum())
for datLab,errLab,label,c in zip(datLabs,errLabs,labels,colors):
    y1 = df['cool',datLab,'rho90kpc']
    y2 = df['cool',datLab,'rho90Rvir']
    err1 = errors[errLab,'stdpkpc']
    err2 = errors[errLab,'stdRvir']
    
    ax1.plot(z,y1,label=label,color=c)
    ax2.plot(z,y2,label=label,color=c)

    lower = y1-err1
    upper = y1+err1
    ax1.fill_between(z,upper,lower,color=c,alpha=0.5)
       
    lower = y2-err2
    upper = y2+err2
    ax2.fill_between(z,lower,upper,color=c,alpha=0.5)

ax1.set_ylabel(r'$r_{90}$ [pkpc]')
ax2.set_ylabel(r'$r_{90}$ [Rvir]')
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
for ax in (ax1,ax2):
    ax.set_xlabel('Redshift')
    ax.invert_xaxis()
    ax.set_xlim([4,1])
    
fig1.savefig(loc+'filamentSizeErrors_pkpc_fill.pdf',bbox_inches='tight',dpi=300)
fig2.savefig(loc+'filamentSizeErrors_rvir_fill.pdf',bbox_inches='tight',dpi=300)
plt.close(fig1)
plt.close(fig2)


# In[29]:

df.info()


# In[47]:

x = np.arange(0,10)
u = x**1.5
l = np.sqrt(x)


# In[61]:

fig,ax = plt.subplots(1,1)
ax.plot(x,x)
ax.plot(x,u)
ax.plot(x,l)
ax.fill_between(x,l,u,color='purple',alpha=0.5)
np.isfinite(x)


# In[ ]:



lower = df['cool',datLab,'rho90Rvir']-errors[errLab,'stdRvir']
upper = df['cool',datLab,'rho90Rvir']+errors[errLab,'stdRvir']
ax2.fill_between(z,lower,upper)

ax1.set_ylabel(r'$r_{90}$ [pkpc]')
ax2.set_ylabel(r'$r_{90}$ [Rvir]')

for ax in (ax1,ax2):
ax.legend(loc='best')
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_xlim([4,1])

fig1.savefig(loc+'filamentSizeErrors_pkpc_fill.pdf',bbox_inches='tight',dpi=300)
fig2.savefig(loc+'filamentSizeErrors_rvir_fill.pdf',bbox_inches='tight',dpi=300)
plt.close(fig1)
plt.close(fig2)


# In[ ]:



