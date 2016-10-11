
# In[11]:

from __future__ import print_function
import matplotlib as mp
mp.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.linalg as sl
import scipy.stats as st
import sys


zmin = 0
zmax = 6
xmean, ymean = 0, 0
xsigma, ysigma = 2, 2

rayVar, datVar, ratio = [], [], []
rayStd, datStd, ratioS = [], [], []

numpoints = 10**6
x = np.random.normal(loc=xmean, scale=xsigma, size=numpoints)
y = np.random.normal(loc=ymean, scale=ysigma, size=numpoints)
z = np.random.uniform(low=zmin, high=zmax, size=numpoints)

# Make datafram
df = pd.DataFrame({'x':x, 'y':y, 'z':z})

# Fit line of best fit
locM = df - df.mean()
covar = locM.cov()
(u,s,v) = sl.svd(covar)
u0 = u[0,2]
u1 = u[1,2]
u2 = u[2,2]
print(u)
print(s)
print(v)

# Get distance ot this line for each cell
locM['a0'] = u2*locM['y'] - u1*locM['z']
locM['a1'] = u0*locM['z'] - u2*locM['x']
locM['a2'] = u1*locM['x'] - u0*locM['y']
locM['dist'] = np.sqrt(locM['a0']**2 + locM['a1']**2 + locM['a2']**2)
locM['r'] = np.sqrt(locM['x']**2 + locM['y']**2)

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

# Get standard deviation
rs = st.rayleigh.std(loc=loc,scale=scale)
ds = locM['dist'].std()
rayStd.append(rs)
datStd.append(ds)
ratioS.append(rs/ds)

# Try for directly measuring distance
paramdir = st.rayleigh.fit(locM['r'])
loc = paramdir[0]
scale = paramdir[1]
rvdir = st.rayleigh.var(loc=loc,scale=scale)
dvdir = locM['dist'].var()
rsdir = st.rayleigh.std(loc=loc,scale=scale)
dsdir = locM['dist'].std()
print('Fit line: scale = {0:.3f}\tvar = {1:.3f}\tstd = {2:.3f}'.format(param[1],rv,rs))
print('Directly: scale = {0:.3f}\tvar = {1:.3f}\tstd = {2:.3f}'.format(paramdir[1],rvdir,rsdir))
print('Compare to data (fit):\tvar = {0:.3f}\tstd = {1:.3f}'.format(locM['dist'].var(),locM['dist'].std()))
print('Compare to data (dir):\tvar = {0:.3f}\tstd = {1:.3f}'.format(locM['r'].var(),locM['r'].std()))

lenu = np.sqrt(u0**2 + u1**2 + u2**2)
angle = np.degrees(np.arccos(u2/len(u)))
print(lenu)
print(angle)

fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
ax1.scatter(x,y,marker='o',alpha=0.01)
ax2.scatter(x,z,marker='o',alpha=0.01)
ax3.scatter(y,z,marker='o',alpha=0.01)

t = np.linspace(-1000,1000,1000)
xp,yp,zp = [],[],[]
for i in range(len(t)):
    xp.append(locM['x'].mean() + t[i]*u0)
    yp.append(locM['y'].mean() + t[i]*u1)
    zp.append(locM['z'].mean() + t[i]*u2)

ax1.plot(xp,yp,'r')
ax2.plot(xp,zp,'r')
ax3.plot(yp,zp,'r')

ax1.set_xlim([x.min(),x.max()])
ax2.set_xlim([x.min(),x.max()])
ax3.set_xlim([y.min(),y.max()])
ax1.set_ylim([y.min(),y.max()])
ax2.set_ylim([z.min(),z.max()])
ax3.set_ylim([z.min(),z.max()])

ax1.set_xlabel('x')
ax2.set_xlabel('x')
ax3.set_xlabel('y')
ax1.set_ylabel('y')
ax2.set_ylabel('z')
ax3.set_ylabel('z')

fig.tight_layout()
fig.savefig('model.png',bbox_inches='tight')



sys.exit()



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

