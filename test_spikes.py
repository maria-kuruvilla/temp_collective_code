# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 11:31:33 2020

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
trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/9/4/session_GS_4_T_9_2/trajectories_wo_gaps/trajectories_wo_gaps.npy'
tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length')
tr.new_time_unit(tr.params['frame_rate'], 'seconds')

frame_range = range(tr.s.shape[0])


for focal in range(tr.number_of_individuals):
    plt.plot(np.asarray(frame_range),tr.speed[frame_range, focal])
    #plt.plot(np.asarray(frame_range),tr.acceleration[frame_range, focal])


tr_dict = np.load(trajectories_file_path, allow_pickle=True).item()
estimated_accuracy = np.nanmean(tr_dict['id_probabilities'])
print(estimated_accuracy)

c = [i for i, value in enumerate(tr.speed[:,2]) if value > 5]

def spikes(trajectory):
    list1 = []
    for j in range(trajectory.number_of_individuals):
        list1 = list1 + [i for i, value in enumerate(trajectory.speed[:,j]) if value > 5]
    return(len(list1)/trajectory.number_of_individuals)

def spikes_position(trajectory):
    list1 = []
    for j in range(trajectory.number_of_individuals):
        list1 = list1 + [i for i, value in enumerate(trajectory.speed[:,j]) if value > 5]
    return(list1)
        
rows = []
with open('looms.csv', 'r') as csvfile:
    looms = csv.reader(csvfile)
    for row in looms:
        rows.append(row)
        
 [(ix,iy) for ix, r in enumerate(rows) for iy, i in enumerate(r) if i == '16']
    
        
            
    
