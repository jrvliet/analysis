
'''
Plots all spectra for a given LOS
Takes in galID< expn, and inclination from command line
'''

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys

# Check command line arguements
if len(sys.argv)!=4:
    print('Usage: python allSpectra.py <galNum> <expn> <inclination>')
    print('Example: python allSpectra.py 27 0.400 90')
    print('Only works with vela2b simulations')
    sys.exit()
else:
    galNum = sys.argv[1]
    expn = sys.argv[2]
    inc = sys.argv[3]

baseLoc = '/mnt/cluster/abs/cgm/vela2b/'
subLoc = 'vela{0:d}/a{1:s}/i{2:s}/'

loc = baseLoc+subLoc.format(galNum,expn,inc)

# Get the ions
with open(loc+'mockspec.config') as f:
    config = f.readlines()

ionSection = config[26:]
ions = []
for line in ionSection:
    ion = line.split()[0]
    ions.append(ion)

print(ions)





