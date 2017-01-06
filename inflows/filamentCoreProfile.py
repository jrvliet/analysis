
# coding: utf-8

# In[34]:

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import statsmodels.api as sm
import sys


# In[17]:

loc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'


# In[12]:

fname = '/home/jacob/research/code/analysis/inflows/filamentMassSize_simple.h5'
dist = pd.read_hdf(fname,'data')

expns = np.arange(0.200,0.500,0.01)


num = 50
inter,slope = [],[]
for a in expns:
    fname = loc+filename.format(a)
    df = pd.read_hdf(fname,'data')
    aname = 'a{0:d}'.format(int(a*1000))
    edge = dist['cool','dense','rho90kpc'][aname]
    df['rho'] = np.sqrt(df['xRot']**2 + df['yRot']**2)
    spaceInd = (df['rho']<edge) & (df['theta']<80)
    denseInd = (df['density']<10**-2.5) & (df['density']>10**-5.5)
    tempInd = (df['temperature']<10**4.5) & (df['temperature']>10**3.5)
    index = spaceInd & denseInd & tempInd
    fil = df[index]
    
    rhoBins = np.linspace(0,edge,num+1)
    r,n = np.zeros(num),np.zeros(num)
    for i in range(num):
        rhoMin = rhoBins[i]
        rhoMax = rhoBins[i+1]
        index = (fil['rho']>rhoMin) & (fil['rho']<rhoMax)
        f = fil[index]
        r[i] = (rhoMin+rhoMax)/2.
        n[i] = np.log10(f['density'].mean())
    x = sm.add_constant(r)
    model = sm.OLS(n,x)
    results = model.fit()
    inter.append(results.params[0])
    slope.append(results.params[1])
    print(a,edge)

    fig,ax = plt.subplots(1,1,figsize=(5,5))
    ax.plot(r,n)
    y = results.params[0] + results.params[1]*r
    print(len(y),len(r))
    ax.plot(r,y)
    ax.set_xlabel('Rho')
    ax.set_ylabel('Density')
    ax.set_title(a)
    s = 'filamentCoreProfile_a{0:.3f}.png'.format(a)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)


    


fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.plot(expns,inter)
ax1.set_xlabel('Expn')
ax1.set_ylabel('Intercept')

ax2.plot(expns,slope)
ax2.set_xlabel('Expn')
ax2.set_ylabel('Slope')

fig.tight_layout()
s = 'filamentCore_fitParams.png'
fig.savefig(s,bbox_inches='tight',dpi=300)

sys.exit()




rhoBins = np.linspace(0,edge,num+1)
r,n = np.zeros(num),np.zeros(num)
for i in range(num):
    rhoMin = rhoBins[i]
    rhoMax = rhoBins[i+1]
    index = (fil['rho']>rhoMin) & (fil['rho']<rhoMax)
    f = fil[index]
    r[i] = (rhoMin+rhoMax)/2.
    n[i] = np.log10(f['density'].mean())


# In[54]:

x = sm.add_constant(r)
model = sm.OLS(n,x)
results = model.fit()


# In[55]:

results.params


# In[56]:

