
'''
Plots the results from vradFraction.py
'''

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt


dataloc = '/home/jacob/research/code/analysis/inflows/'

vrMean = pd.read_hdf(dataloc+'vradFraction_mean.h5','data')
vrStd = pd.read_hdf(dataloc+'vradFraction_std.h5','data')


rMean = pd.read_hdf(dataloc+'r_mean.h5','data')
rStd = pd.read_hdf(dataloc+'r_std.h5','data')


fig,axes = plt.subplots(2,2,figsize=(10,10))
axes[0].plot(vrMean


