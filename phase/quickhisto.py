
import numpy as np
import matplotlib.pyplot as plt


location = '/home/matrix3/jrvander/sebass_gals/dwarfs/abscells/'
filename = 'D9m4a.CIV.bulk_abscells_extend.dat'

data = np.loadtxt(location+filename, skiprows=1)

temp = data[:,8]
nH = data[:,7]

plt.hist(temp, 50, histtype='step')
plt.title('Temperature')
plt.savefig('D9m4a.CIV.temp.hist.pdf')

plt.cla()
plt.clf()

plt.hist(nH, 50, histtype='step')
plt.title('Density')
plt.savefig('D9m4a.CIV.dense.hist.pdf')

