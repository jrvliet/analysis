
from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt

filename = 'vela2b_cloud_loc_stats.h5'
df = pd.read_hdf(filename, 'data')

fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(15,5))
axes = (ax1,ax2,ax3)
ax1a = ax1.twinx()
ax2a = ax2.twinx()
ax3a = ax3.twinx()

lns1 = ax1.plot(df['galNum'], df['r1'], 'xr', label='Radius')
lns2 = ax1a.plot(df['galNum'], df['theta1'], 'og', label='Theta')
lns3 = ax1a.plot(df['galNum'], df['phi1'], '^b', label='Phi')
ax1.set_xlabel('Gal Number')
ax1.set_ylabel('Average Property [kpc]')
ax1a.set_ylabel('Average Property [deg]')
ax1.set_title('Cloud 1')
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, frameon=False)

lns1 = ax2.plot(df['galNum'], df['r2'], 'xr', label='Radius')
lns2 = ax2a.plot(df['galNum'], df['theta2'], 'og', label='Theta')
lns3 = ax2a.plot(df['galNum'], df['phi2'], '^b', label='Phi')
ax2.set_xlabel('Gal Number')
ax2.set_ylabel('Average Property [kpc]')
ax2a.set_ylabel('Average Property [deg]')
ax2.set_title('Cloud 2')
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs, frameon=False)

lns1 = ax3.plot(df['galNum'], df['r3'], 'xr', label='Radius')
lns2 = ax3a.plot(df['galNum'], df['theta3'], 'og', label='Theta')
lns3 = ax3a.plot(df['galNum'], df['phi3'], '^b', label='Phi')
ax3.set_xlabel('Gal Number')
ax3.set_ylabel('Average Property [kpc]')
ax3a.set_ylabel('Average Property [deg]')
ax3.set_title('Cloud 3')
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax3.legend(lns, labs, frameon=False)



fig.tight_layout()
fig.savefig('vela2b_cloud_averageProps.png', bbox_inches='tight', dpi=300)

    


