# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 13:09:50 2020

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

sigma_values = 1.5
trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/9/16/session_GS_16_T_9_2/trajectories_wo_gaps/trajectories_wo_gaps.npy'
#trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/29/1/session_GS_1_T_29_2/trajectories/trajectories.npy'
trajectory = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length')

nnd = np.empty([trajectory.number_of_individuals, 1])
nnd.fill(np.nan)
nd = np.empty([trajectory.number_of_individuals, 1])
nd.fill(np.nan)
annd_per_frame = np.zeros([trajectory.s.shape[0], 1])
#annd_per_frame = np.empty([5, 1])
annd_per_frame.fill(np.nan)
k = 0
for frame in range(trajectory.s.shape[0]):
        
        for i in range(trajectory.number_of_individuals):
            
            for j in range(trajectory.number_of_individuals):
                
                if i!=j:
                    nd[k] = distance.euclidean(trajectory.s[frame,i,:], trajectory.s[frame,j,:]) # distance between every individual and every other individual (neighbour distance)
                    k = k+1
                    
                    
            nnd[i] = np.nanmin(nd) #nearest neighbor distance
            
            k = 0
        annd_per_frame[frame] = np.nanmean(nnd) #average nearest neighbour distance per fram
        
annd = np.nanmean(annd_per_frame) 
        
                
print(annd)
