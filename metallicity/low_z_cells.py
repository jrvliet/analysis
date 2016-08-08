
import pandas as pd


cloudTlow = 10**4
cloudThi = 10**4.5
cloudNlow = 10**-5
cloudNhi = 10**-4

with open('galaxy.props', 'r') as f:

    galID = f.readline().split()[1]
    expn = f.readline().split()[1]
    redshift = f.readline().split()[1]
    mvir = float(f.readline().split()[1])
    rvir = float(f.readline().split()[1])
    inc = int(float(f.readline().split()[1]))


boxname = '{0:s}_GZa{1:s}.h5'.format(galID,expn)
box = pd.read_hdf(boxname, 'data')

snIICut = box.quantile(0.1)['SNII']
print snIICut

lowZname = '{0:s}_a{1:s}_losZcells.txt'.format(galID, expn)
cloudname = '{0:s}_a{1:s}_cloudcells.txt'.format(galID, expn)
flowZ = open(lowZname, 'w')
fcloud = open(cloudname, 'w')

for i in range(len(box['SNII'])):
    
    if box['SNII'][i]<=snIICut:
        flowZ.write('{0:d}\n'.format(i+1))
    
    if (box['density'][i]<cloudNhi and box['density'][i]>cloudNlow and 
        box['temperature'][i]>cloudTlow and box['temperature'][i]<cloudThi):
        fcloud.write('{0:d}\n'.format(i+1))

flowZ.close()
fcloud.close()


