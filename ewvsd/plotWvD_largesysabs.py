import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages

allfile = sys.argv[1]

# Parse the allfile
galID = allfile.split('.')[0]
if 'master' in allfile:
    a = 'master'
    ion = allfile.split('.')[1]
else:
    a1 = allfile.split('.')[2].split('a')[1]
    a2 = allfile.split('.')[3]
    a = a1+'.'+a2
    ion = allfile.split('.')[1]

# Read in the ALL file
ew_raw, imp_raw  = np.loadtxt(allfile, skiprows=1, usecols=(5, 22), unpack=True)

# Remove any EW that is zero, as these come from LOS with
# no detections
ew = []
imp = []
for i in range(0,len(ew_raw)):
    if ew_raw[i]!=0.0:
        ew.append(ew_raw[i])
        imp.append(imp_raw[i])

ew = np.array(ew)
imp = np.array(imp)
ew = np.log10(ew)

Dmin =  []
Dmax = []
for i in range(0,15):
    Dmin.append(i*0.1)
    Dmax.append((i+1)*0.1)
    
ew_mean = []
ew_err = []
imp_mean = []
imp_err_pos = []
imp_err_neg = []
for i in range(0,len(Dmin)):
    
    low = Dmin[i]
    high = Dmax[i]
    
    ewtotal = []
    imptotal = []
    for j in range(0,len(imp)):
        
        if imp[j]<high and imp[j]>low:
            ewtotal.append(ew[j])
            imptotal.append(imp[j])
            
    im = np.mean(imptotal)
    imp_mean.append(im)
    ew_mean.append(np.mean(ewtotal))
    ew_err.append(np.std(ewtotal))
    imp_err_pos.append(high-im)
    imp_err_neg.append(im-low)
    

plt.errorbar(imp_mean, ew_mean, yerr=ew_err, xerr=[imp_err_neg, imp_err_pos], color='black', fmt='s', linewidth=2)
plt.savefig('ew_d_quick_binned_'+galID+'_'+a+'_'+ion+'.pdf')
