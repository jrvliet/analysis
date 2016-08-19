
import pandas as pd

import matplotlib.pyplot as plt







haloloc = '/mnt/cluster/abs/Simulations/vela2.1/VELA{0:s}/output/ana/'
haloname = 'halos_{0:s}.txt'

gasloc = '/mnt/cluster/abs/cgm/vela2b/vela{0:d}/{1:s}/'
gasname = 'vela2b-{0:d}_GZa{1:s}.h5'

galNums = range(21,30)
expns = ['0.490']*len(galNums)
expns[galNums.index(24)] = '0.450'

for galNum, a in zip(galNums,expns):


    hloc = haloloc.format(galNum)
    hfile = haloname.format(a)
    gloc = gasloc.format(galNum,a)
    gfile = gasname.format(galNum,a)


    df = pd.read_hdf(gloc+gfile, 'data')

    hx, hy, hz, hrvir = np.loadtxt(hloc+hfile, skiprows=2, usecols=(15,16,17,8)

    




