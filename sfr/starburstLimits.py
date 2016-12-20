
# coding: utf-8

# In[1]:

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[2]:

loc = '/home/jacob/research/code/analysis/sfr/'
filename = 'vela2b-27_sfr.csv'
fname = loc+filename


# In[3]:

df = pd.read_csv(fname)


# In[13]:

sigma = df['sfr'].std()
ave = df['sfr'].mean()
med = df['sfr'].median()


# In[14]:

med


# In[7]:

sigma


# In[10]:

ave


# In[24]:

baseSFR = df['sfr'][df['sfr']<med+sigma].mean()
baseSFRStd = df['sfr'][df['sfr']<med+sigma].std()


# In[21]:

df['sfr'][df['sfr']>baseSFR+sigma]


# In[22]:

df['sfr'][df['sfr']>med+sigma]


# In[23]:

df['sfr'][df['sfr'] >ave+sigma]


# In[25]:

baseSFRStd


# In[33]:

print(df[['a','sfr']][df['sfr']>baseSFR+baseSFRStd])
print(df[['a','sfr']][df['sfr']>med+baseSFRStd])
print(df[['a','sfr']][df['sfr']>ave+baseSFRStd])


# In[34]:

burst = df[['a','sfr']][df['sfr']>ave+baseSFRStd]


# # Define outliers using quantiles

# In[42]:

q1 = np.percentile(df['sfr'],25)
q3 = np.percentile(df['sfr'],75)
iqr = q3-q1


# In[43]:

lq,uq,iqr


# In[44]:

lowerLimit = q1 - 1.5*iqr
upperLimit = q3 + 1.5*iqr


# In[50]:

baseInds = df['sfr']<upperLimit


# In[49]:

df[['a','sfr']][df['sfr']<upperLimit].mean()


# ### Get the baseline SFR without the outliers 

# In[53]:

baseSFRMean = df['sfr'][baseInds].mean()
baseSFRStd = df['sfr'][baseInds].std()


# ### Select out the burst 
# Selecting points with SFR greater than a multiple of the standard deviation from the mean, using only the baseline mean and standard deviations (ie, no outliers)

# In[ ]:

multiple = 1
burst = df[['a','sfr']][df['sfr']>baseSFRMean+multiple*baseSFRStd]


# In[56]:

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(df['a'],df['sfr'])
ax.plot(burst['a'],burst['sfr'],linestyle='none',marker='s',color='r')


# In[57]:

burst['a']


# In[ ]:



