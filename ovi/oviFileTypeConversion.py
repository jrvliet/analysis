import pandas as pd

loc = '/lustre/projects/p089_swin/jvander/analysis/ovi/'
filename = 'vela2b-{0:d}.a0.490.i90_OVIcellsGalFrame.h5'
galNums = range(21,30)

for galNum in galNums:
    fname = loc+filename.format(galNum)
    try:
        df = pd.read_hdf(fname,'data')
    except IOError:
        continue
    outname = fname.replace('h5','csv')
    df.to_csv(outname,index=False)
