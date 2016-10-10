
# In[11]:

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.linalg as sl
import scipy.stats as st




zmin = 0
zmax = 3
xmean, ymean = 0, 0
xsigma, ysigma = 2, 2


# In[25]:

rayVar, datVar, ratio = [], [], []
rayStd, datStd, ratioS = [], [], []

nums = range(1,9)
for num in nums:
    print(num)
    numpoints = 10**num
    x = np.random.normal(loc=xmean, scale=xsigma, size=numpoints)
    y = np.random.normal(loc=ymean, scale=ysigma, size=numpoints)
    z = np.random.uniform(low=zmin, high=zmax, size=numpoints)


    # Make datafram
    df = pd.DataFrame({'x':x, 'y':y, 'z':z})

    # Fit line of best fit
    locM = df - df.mean()
    covar = locM.cov()
    (u,s,v) = sl.svd(covar)
    u0 = u[0,0]
    u1 = u[1,0]
    u2 = u[2,0]

    # Get distance ot this line for each cell
    locM['a0'] = u2*locM['y'] - u1*locM['z']
    locM['a1'] = u0*locM['z'] - u2*locM['x']
    locM['a2'] = u1*locM['x'] - u0*locM['y']
    locM['dist'] = np.sqrt(locM['a0']**2 + locM['a1']**2 + locM['a2']**2)

    # Fit rayleigh to distance distribution
    param = st.rayleigh.fit(locM['dist'])
    loc = param[0]
    scale = param[1]


    # Get variance
    rv = st.rayleigh.var(loc=loc,scale=scale)
    dv = locM['dist'].var()
    rayVar.append(rv)
    datVar.append(dv)
    ratio.append(rv/dv)
    print(rv/dv)
    # Get standard deviation
    rs = st.rayleigh.std(loc=loc,scale=scale)
    ds = locM['dist'].std()
    rayStd.append(rs)
    datStd.append(ds)
    ratioS.append(rs/ds)
    print(rv/dv)
    #xfit = np.linspace(locM['dist'].min(),locM['dist'].max())
    #yfit = st.rayleigh.pdf(xfit, loc=loc, scale=scale)
    #fig,ax = plt.subplots(1,1,figsize=(5,5))
    #ax.hist(locM['dist'],bins=500,histtype='step',normed=True)
    #ax.plot(xfit,yfit)



fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))
ax1.plot(nums,rayVar,label='Fit')
ax1.plot(nums,datVar,label='Data')
ax1.set_xlabel('log Number of Points')
ax1.set_ylabel('Variance')

ax2.plot(nums,ratio)
ax2.set_xlabel('log Number of Points')
ax2.set_ylabel('Fit Variance / Data Variance')
ax2.set_ylim([0,1])

ax3.plot(nums,rayStd,label='Fit')
ax3.plot(nums,datStd,label='Data')
ax3.set_xlabel('log Number of Points')
ax3.set_ylabel('Standard Deviation')

ax4.plot(nums,ratioS)
ax4.set_xlabel('log Number of Points')
ax4.set_ylabel('Fit Std Dev / Data Std Dev')
ax4.set_ylim([0,1])

fig.tight_layout()
s = 'model_inflow_variance.png'
fig.savefig(s,bbox_inches='tight',dpi=300)

