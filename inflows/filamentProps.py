
# coding: utf-8

# In[1]:

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')

pd.options.mode.chained_assignment = None

def mfToZ(Z_cell):

    # mfToZ.py
    # 
    # Converts mass fractions to Z/Z_sun
    #
    # Z/Z_sun = (Z_m / X_h)_cell
    #           ----------------
    #           (Z_m / X_h)_sun
    #
    # where X_h + Y_he + Z_m = 1
    # Since the simulation does not track Helium, need to assume 
    # a value for r = Y/X
    

    # Solar values
    # Taken from Chris' Notes 
    X_sun = 0.70683
    Y_sun = 0.27431
    Z_sun = 0.0188
    
    r = 0.3347

    # Loop through cells
    X_cell = (1 - Z_cell) / (1 + r)
    Z = (Z_cell / X_cell) / (Z_sun / X_sun)

    return Z

# In[ ]:

def mkLines(ax):
    points = [2.57,2.13,1.86]
    ymin,ymax = ax.get_ylim()
    for point in points:
        ax.vlines(point,ymin,ymax,linestyle='dashed')
    ax.set_ylim([ymin,ymax])


# In[2]:

dataloc = '/home/jacob/research/velas/vela2b/vela27/'
filename = 'a{0:.3f}/vela2b-27_GZa{0:.3f}.rot.h5'
expns = np.arange(0.200,0.550,0.01)
expns = np.arange(0.200,0.500,0.01)


# In[3]:

loN,hiN = -5.5,-2.5
loT,hiT = 3.5,4.5


# In[4]:

kpc2km = 3.086e16 
pc2cm = 3.086e18 
mH = 1.6737e-24 
mSun = 1.989e33 
s2yr = 3.154e7


readIn = 1
# In[5]:
hdfFile = 'filamentProps.h5'
if readIn == 0:
    mass,dense,temp,redshift,metal,speed = [],[],[],[],[],[]
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
        cloud['metal'] = mfToZ(cloud['SNII']+cloud['SNIa'])
        cloud['speed'] = np.sqrt(cloud['vx']**2 + cloud['vy']**2 +
                                 cloud['vz']**2)

        totalMass = cloud['mass'].sum()
        meanDense = cloud['density'].median()
        medianTemp = cloud['temperature'].median()
        medianSpeed = cloud['speed'].median()
        medianMetal = cloud['metal'].median()

        redshift.append(1./a - 1)
        mass.append(np.log10(totalMass))
        dense.append(np.log10(meanDense))
        temp.append(np.log10(medianTemp))
        metal.append(np.log10(medianMetal))
        speed.append(medianSpeed)
        

    df = pd.DataFrame()
    df['totalMass'] = mass
    df['medianDense'] = dense
    df['medianTemp'] = temp
    df['redshift'] = redshift
    df['expn'] = expns
    df['speed'] = speed
    df['metal'] = metal
    df.to_hdf(hdfFile,'data')

else:
    df = pd.read_hdf(hdfFile,'data')
    redshift = df['redshift']
    mass = df['totalMass']
    dense = df['medianDense']
    temp = df['medianTemp']
    speed = df['speed']
    metal = df['metal']


# In[10]:

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(redshift,mass)
mkLines(ax)
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_ylabel('Total Mass [log Msun]')
fig.savefig('filamentProps_mass.pdf',bbox_inches='tight',dpi=300)


# In[ ]:

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(redshift,dense)
ax.set_ylim([-5.5,-2.5])

# Plot lines to demark the different density bins
cuts = [-3.5,-4.5]
xmin,xmax = ax.get_xlim()
for cut in cuts:
    ax.hlines(cut,xmin,xmax,linestyle='--',color='g')
x = np.linspace(xmin,xmax,1000)
alpha = 0.15
ax.fill_between(x,-4.5,-5.5,color='red',alpha=alpha)
ax.fill_between(x,-3.5,-4.5,color='green',alpha=alpha)
ax.fill_between(x,-2.5,-3.5,color='blue',alpha=alpha)

ax.text(3.8,-4.7,r'\textit{diffuse}',color='red')
ax.text(3.8,-3.7,r'\textit{mid}',color='green')
ax.text(3.8,-2.7,r'\textit{dense}',color='blue')

mkLines(ax)
ax.set_xlim([xmin,xmax])
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_ylabel('Median Density [log cm$^{-3}$]')
fig.savefig('filamentProps_density.pdf',bbox_inches='tight',dpi=300)
plt.close(fig)

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(redshift,temp)
mkLines(ax)
ax.set_ylim([3.5,4.5])
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_ylabel('Median Temperature [log K]')
fig.savefig('filamentProps_temperature.pdf',bbox_inches='tight',dpi=300)
plt.close(fig)

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(redshift,speed)
mkLines(ax)
#ax.set_ylim([3.5,4.5])
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_ylabel('Median Speed [km/s]')
fig.savefig('filamentProps_speed.pdf',bbox_inches='tight',dpi=300)
plt.close(fig)

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(redshift,metal)
mkLines(ax)
#ax.set_ylim([3.5,4.5])
ax.set_xlabel('Redshift')
ax.invert_xaxis()
ax.set_ylabel('Median Metallicity')
fig.savefig('filamentProps_metallicity.pdf',bbox_inches='tight',dpi=300)
plt.close(fig)

