# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:42:34 2020

@author: Maria Kuruvilla
"""

import os
import pathlib
from pprint import pprint

import numpy as np
from scipy import stats
from scipy.spatial import distance
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import trajectorytools as tt
import trajectorytools.plot as ttplot
import trajectorytools.socialcontext as ttsocial
from trajectorytools.constants import dir_of_data
import pandas
import csv

sigma_values = 1.5
#trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/29/4/session_GS_4_T_29_1/trajectories/trajectories.npy'

trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/9/4/session_GS_4_T_9_2/trajectories_wo_gaps/trajectories_wo_gaps.npy'
tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length')
tr.new_time_unit(tr.params['frame_rate'], 'seconds')

frame_range = range(tr.s.shape[0])


for focal in range(tr.number_of_individuals):
    plt.plot(np.asarray(frame_range),tr.speed[frame_range, focal])
    #plt.plot(np.asarray(frame_range),tr.acceleration[frame_range, focal])

plt.xlabel('Frame number')
plt.ylabel('Speed (BL/s)')


rows = []
with open('looms.csv', 'r') as csvfile:
    looms = csv.reader(csvfile)
    for row in looms:
        rows.append(row)
        
loom = np.zeros([5,1])        
for i in range(len(rows)):
    if rows[i][1]=='Cam 12' and rows[i][3]=='1' and rows[i][4]=='1':
        for j in range(5):
            loom[j] = int(rows[i][2]) + j*11403
        
        
        
        