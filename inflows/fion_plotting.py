
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

loc = '/lustre/projects/p089_swin/jvander/analysis/inflows/'
filename = 'fIon.h5'
fname = loc+filename

df = pd.read_hdf(fname,'data')




