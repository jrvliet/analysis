
# coding: utf-8

# In[1]:

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')


# In[ ]:

def mkLines(ax):
    ymin,ymax = ax.get_ylim()
    ax.vlines(0.29,ymin,ymax,linestyle='dashed')
    ax.vlines(0.32,ymin,ymax,linestyle='dashed')
    ax.vlines(0.35,ymin,ymax,linestyle='dashed')
    ax.set_ylim([ymin,ymax])


# In[2]:

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
expns = np.arange(0.200,0.550,0.01)


# In[3]:

loN,hiN = -5.5,-2.5
loT,hiT = 3.5,4.5


# In[4]:

kpc2km = 3.086e16 
pc2cm = 3.086e18 
mH = 1.6737e-24 
mSun = 1.989e33 
s2yr = 3.154e7


# In[5]:

mass,dense,temp,redshift = [],[],[],[]
for a in expns:
    print(a)
    fname = dataloc+filename.format(a)
    df = pd.read_hdf(fname,'data')
    
    tempInd = (df['temperature']>10**loT) & (df['temperature']<10**hiT)
    densInd = (df['density']>10**loN) & (df['density']<10**hiN)
    spacInd = (df['theta']<80)
    
    cloud = df[tempInd & densInd & spacInd]
    #print(a,len(cloud))
    cloud['mass'] = cloud['density']*mH*(cloud['cell_size']*pc2cm)**3 / mSun

    totalMass = cloud['mass'].sum()
    meanDense = cloud['density'].median()
    medianTemp = cloud['temperature'].median()
    mass.append(np.log10(totalMass))
    dense.append(np.log10(meanDense))
    temp.append(np.log10(medianTemp))
    redshift.append(1./a - 1)


# In[10]:

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(redshift,mass)
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_ylabel('Total Mass [log Msun]')
fig.savefig('filamentProps_mass.pdf',bbox_inches='tight',dpi=300)


# In[ ]:

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(redshift,dense)
ax.set_ylim([-5.5,-2.5])
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_ylabel('Median Density [log cm$^{-3}$]')
fig.savefig('filamentProps_density.pdf',bbox_inches='tight',dpi=300)


# In[ ]:

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(redshift,temp)
ax.set_ylim([3.5,4.5])
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_ylabel('Median Temperature [log K]')
fig.savefig('filamentProps_temperature.pdf',bbox_inches='tight',dpi=300)

