
'''
Plots all spectra for a given LOS
Takes in galID< expn, and inclination from command line
'''

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob

# Check command line arguements
if len(sys.argv)!=5:
    print('Usage: python allSpectra.py <galNum> <expn> <inclination> <losnum>')
    print('Example: python allSpectra.py 27 0.400 90 52')
    print('Only works with vela2b simulations')
    sys.exit()
else:
    galNum = sys.argv[1]
    expn = sys.argv[2]
    inc = sys.argv[3]
    los = int(sys.argv[4])

baseLoc = '/mnt/cluster/abs/cgm/vela2b/'
subLoc = 'vela{0:s}/a{1:s}/i{2:s}/'

loc = baseLoc+subLoc.format(galNum,expn,inc)

# Get the ions
with open(loc+'mockspec.config') as f:
    config = f.readlines()

ionSection = config[26:]
ions = []
for line in ionSection:
    ion = line.split()[0]
    ions.append(ion)


for ion in ions:
    query = '{0:s}/{1:s}/vela2b-{2:s}.{1:s}.los{3:04d}.*.spec'.format(loc,ion,galNum,los)
    specs = glob.glob(query)

    print(ion)
    print(specs)
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    for specfile in specs:
        transition = specfile.split('.')[-2]
        data = np.loadtxt(specfile)
        wavelength = data[:,0]
        velocity = data[:,1]
        flux = data[:,2]
        ax.plot(velocity,flux,label=transition)
    ax.legend(loc='best')
    ax.set_ylim([0,1.1])
    ax.set_title(ion)
    s = 'vela2b-{0:s}_a{1:s}_i{2:s}_los{3:04d}_{4:s}.png'.format(galNum,expn,inc,los,ion)
    fig.savefig(s,bbox_inches='tight',dpi=300)
    plt.close(fig)
    




