
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st

fig, ax = plt.subplots(figsize=(5,5))
x = np.linspace(0,100,10000)

loc = 0
scales = range(1,10)

for scale in scales:
    
    y = st.rayleigh.pdf(x, loc=loc, scale=scale)
    ax.plot(x,y,label='{0:d}'.format(scale))

ax.legend()
fig.savefig('rayleigh_scales.png',bbox_inches='tight',dpi=300)


