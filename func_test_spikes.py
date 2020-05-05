# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 10:56:39 2020

@author: Maria Kuruvilla
"""


#function to input trajectory and get speed with GS,T,rep.


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
trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/9/4/session_GS_4_T_9_2/trajectories_wo_gaps/trajectories_wo_gaps.npy'
tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length')
tr.new_time_unit(tr.params['frame_rate'], 'seconds')

frame_range = range(tr.s.shape[0])

def track_check(tr, temp, group, rep):
    frame_range = range(tr.s.shape[0])
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for focal in range(tr.number_of_individuals):
        ax.plot(np.asarray(frame_range),tr.speed[frame_range, focal])
        
    
    ax.set_xlabel('Frame number')
    ax.set_ylabel('Speed (BL/s)')
    ax.set_title('Temp:' + str(temp) + ' Group:' + str(group) + ' Replicate:' + str(rep))
    return(ax)

temperature = range(9,30,4)



group = [1,2,4,8,16]



replication = range(3)

ii = 0

for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        for k in replication:
            
            
            if j==1:
                try:
                    trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k+1)+'/trajectories/trajectories.npy'
                    
                except FileNotFoundError:
                    print('File not found')
                    continue
            
                    
                
            else:
                try:
                    trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k+1)+'/trajectories_wo_gaps/trajectories_wo_gaps.npy'
            
                except FileNotFoundError:
                    print('File not found')
                    continue
                
            
            sigma_values = 1.5 #smoothing parameter
            tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length') # normalizing by body length
            tr.new_time_unit(tr.params['frame_rate'], 'seconds') # changing time unit to seconds
            track_check(tr, i, j, k+1)
            
        jj= jj + 1
        
    ii = ii + 1