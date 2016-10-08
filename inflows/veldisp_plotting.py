
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[2]:

fname = '/home/jacob/research/code/analysis/inflows/velocity_dispersions.h5'
with pd.HDFStore(fname) as store:
    print store.keys()


# In[4]:

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(15,10))
store = pd.HDFStore(fname)
for key in store.keys():
    df = store.select(key)
    df.plot(x='z',y='sigma_valong',ax=ax1,label=key)
    df.plot(x='z',y='sigma_vperp',ax=ax2,label=key)
    df.plot(x='z',y='sigma_vr',ax=ax4,label=key)
    df.plot(x='z',y='sigma_vtheta',ax=ax5,label=key)
    df.plot(x='z',y='sigma_vphi',ax=ax6,label=key)
for ax in [ax1,ax2,ax4,ax5,ax6]:
    ax.set_xlim([1,4])
    ax.invert_xaxis()
    ax.legend(loc='upper left', fontsize='small')
ax1.set_ylabel('Valong')
ax2.set_ylabel('Vperp')
ax4.set_ylabel('Vr')
ax5.set_ylabel('Vtheta')
ax6.set_ylabel('Vphi')

fig.tight_layout()
s = 'velocity_dispersions_evolution.png'

fig.savefig(s,bbox_inches='tight',dpi=300)

store.close()
