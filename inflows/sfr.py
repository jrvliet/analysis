

from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dataloc = '/mnt/cluster/abs/Simulations/vela2.1/VELA27/output/ana/'
filename = 'halos_{0:.3f}.txt'

times = pd.read_csv('lookback_time.csv')

print(times)

print(times['age [Gyr]'])

for i in range(len(times)):
    
    



