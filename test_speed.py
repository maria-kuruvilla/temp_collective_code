# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 20:41:18 2020

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
from datetime import datetime

import time




sigma_values = 1.5
trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/9/2/session_GS_2_T_9_2/trajectories_wo_gaps/trajectories_wo_gaps.npy'
#trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/29/1/session_GS_1_T_29_2/trajectories/trajectories.npy'
trajectory = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length')

#nnd = np.empty([trajectory.number_of_individuals, 1])
#nnd.fill(np.nan)
#nd = np.empty([trajectory.number_of_individuals, 1])
#nd.fill(np.nan)
nnd = np.zeros([trajectory.s.shape[0], trajectory.number_of_individuals])
#annd_per_frame = np.empty([5, 1])
nnd.fill(np.nan)
k = 0
frames = 1000
start = time.process_time()
for frame in range(frames):
        #print(frame)
        #print(datetime.now())
        
        for i in range(trajectory.number_of_individuals):
            #print(i)
            #print(datetime.now())
            
            for j in range(trajectory.number_of_individuals):
                #print(j)
                #print(datetime.now())
                
                if i!=j:
                    nd[k] = distance.euclidean(trajectory.s[frame,i,:], trajectory.s[frame,j,:]) # distance between every individual and every other individual (neighbour distance)
                    #nd[k]=np.linalg.norm(trajectory.s[frame,i,:] - trajectory.s[frame,j,:])
                    #a_min_b = trajectory.s[frame,i,:] - trajectory.s[frame,j,:]
                    #nd[k] = np.sqrt(np.einsum('ij,ij->j', a_min_b, a_min_b))
                    k = k+1
                    #print(datetime.now())
                    
                    
            nnd[i] = np.nanmin(nd) #nearest neighbor distance
            
            k = 0
        annd_per_frame[frame] = np.nanmean(nnd) #average nearest neighbour distance per fram
        
annd = np.nanmean(annd_per_frame) 
        
                
print(annd)
print(time.process_time() - start)










sigma_values = 1.5
trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/9/4/session_GS_4_T_9_2/trajectories_wo_gaps/trajectories_wo_gaps.npy'
#trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/29/1/session_GS_1_T_29_2/trajectories/trajectories.npy'
trajectory = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length')

#nnd = np.empty([trajectory.number_of_individuals, 1])
#nnd.fill(np.nan)
#nd = np.empty([trajectory.number_of_individuals, 1])
#nd.fill(np.nan)
nnd = np.zeros([trajectory.s.shape[0], trajectory.number_of_individuals])
#annd_per_frame = np.empty([5, 1])
nnd.fill(np.nan)


start = time.process_time()
nd = np.empty([trajectory.s.shape[0],trajectory.number_of_individuals])
nd.fill(np.nan)

for i in range(trajectory.number_of_individuals):
    for j in range(trajectory.number_of_individuals):
        if i!=j:
            nd[:,j] = np.sqrt((trajectory.s[:,i,0] - trajectory.s[:,j,0])**2 + (trajectory.s[:,i,1] - trajectory.s[:,i,1])**2)
            
    nnd[:,i] = np.nanmin(nd,1)
        
annd = np.nanmean(nnd) 

print(annd)
print(time.process_time() - start)