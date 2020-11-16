"""
Goal - to make a heatmap of speed and acceleration of all the trajectories 
"""

import os
import pathlib
from pprint import pprint

import numpy as np
from numpy import *
from scipy import stats
from scipy.spatial import distance
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import trajectorytools as tt
import trajectorytools.plot as ttplot
import trajectorytools.socialcontext as ttsocial
from trajectorytools.constants import dir_of_data
import csv
import pickle

temperature = range(9,30,4)



group = [1,2,4,8,16]



replication = range(10) # number of replicates per treatment

#output parent directory
parent_dir = '../../output/temp_collective/roi'

ii = 0 # to keep count of temperature

#frames = 5000 #number of frames for which annd is calculated
x = [] #append speed values to this 
y = []
for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
        for k in replication:
            
            input_file = out_dir + str(k+1) + '_nosmooth.p'
            
            try:
                tr = pickle.load(open(input_file, 'rb')) # 'rb is for read binary
                print(i,j,k+1)
            except FileNotFoundError:
                print(i,j,k+1)
                print('File not found')
                continue
            for m in range(tr.speed.shape[1]):
                x= np.r_[x,tr.speed[:,m]]
                y= np.r_[y,tr.acceleration[:,m]]
            
                





hist2d = []
vmin = 0
vmax = 0

hist2d.append(np.histogram2d(x[~np.isnan(x)], y[~np.isnan(y)], [20, 20])[0])

# Plot distributions of positions in the arena
figv, ax_hist = plt.subplots(1,figsize=(15,7), sharey=True, sharex=True)
min_x, max_x = np.nanmin(x), np.nanmax(x)
min_y, max_y = np.nanmin(y), np.nanmax(y)
p = ma.log(hist2d[0][:])
vmax = np.max(p.filled(0)) 
ax_hist.imshow(p.filled(0), interpolation = 'none',vmin=vmin, vmax=vmax,origin = 'lower',extent=[min_x, max_x, min_y, max_y])
#ax_hist.set_aspect('equal')
#ticksx = np.arange(min_x, max_x,50)
#ticksy = np.arange(min_y, max_y,500)
##plt.xticks(ticksx)
#plt.yticks(ticksy)
#,extent=[min_y, max_y, min_x, max_x]
ax_hist.set_aspect(0.004)
ax_hist.set_xlabel('Speed')
ax_hist.set_ylabel('Acceleration')
out_dir = parent_dir = '../../output/temp_collective/roi_figures/speed_acc_heatmap_masked_4.png'
figv.savefig(out_dir, dpi = 300)
plt.show()
