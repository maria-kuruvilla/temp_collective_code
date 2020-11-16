# -*- coding: utf-8 -*-
"""
Created on Mon May  4 09:42:34 2020
edited on 10/08/2020

@author: Maria Kuruvilla
goal - one function to return trajectory and one function to plot the tracking results (positions)
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
import csv

def track_check_position(tr, temp, group, rep): #replicates start from 1
    frame_range = range(tr.s.shape[0])
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for i in range(tr.number_of_individuals):
        ax.plot(tr.s[frame_range,i,0], tr.s[frame_range,i,1])
        
    
    ax.set_xlabel('X (BL)')
    ax.set_ylabel('Y (BL)')
    ax.set_title('Trajectories')
    return(ax)

def trajectory_track_check(i, j , k): #takes replicates starting from 1
    if j == 1:
        trajectories_file_path = '../../data/temp_collective/roi/'+str(i)+'/' +str(j)+'/GS_'+str(j)+'_T_'+str(i)+'_roi_'+str(k)+'/trajectories.npy'
    else:
        trajectories_file_path = '../../data/temp_collective/roi/'+str(i)+'/' +str(j)+'/GS_'+str(j)+'_T_'+str(i)+'_roi_'+str(k)+'/trajectories_wo_gaps.npy'
    sigma_values = 1.5 #smoothing parameter
    tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, 		center=True).normalise_by('body_length')#, smooth_params={'sigma': 	sigma_values}).normalise_by('body_length') # normalizing by body length
    tr.new_time_unit(tr.params['frame_rate'], 'seconds') # changing time unit to seconds
    track_check_position(tr,i,j,k)
    return(tr) 
