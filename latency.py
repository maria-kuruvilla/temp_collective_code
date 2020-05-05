# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:56:48 2020

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
import csv

rows = []
with open('looms.csv', 'r') as csvfile:
    looms = csv.reader(csvfile)
    for row in looms:
        rows.append(row)


def trajectory(i, j , k):
    if j == 1:
        trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k)+'/trajectories/trajectories.npy'
    else:
        trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k)+'/trajectories_wo_gaps/trajectories_wo_gaps.npy'
    sigma_values = 1.5 #smoothing parameter
    tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length') # normalizing by body length
    tr.new_time_unit(tr.params['frame_rate'], 'seconds') # changing time unit to seconds
    return(tr) 

def loom_frame(temp, groupsize, rep):
    if temp == 29:
        cam = 'Cam 7'
    elif temp == 25:
        cam = 'Cam 8'
    elif temp == 17:
        cam = 'Cam 9'
    elif temp == 13:
        cam = 'Cam 10'
    elif temp == 21:
        cam = 'Cam 11'
    elif temp == 9:
        cam = 'Cam 12'
    g = str(groupsize)
    r = str(rep)
    loom = np.zeros([5,1])        
    for i in range(len(rows)):
        if rows[i][1]==cam and rows[i][3]==g and rows[i][4]==r:
            for j in range(5):
                if rows[i][2]=='':
                    loom[j] = int(rows[i-1][2]) + 600 + j*11403
                else:
                   loom[j] = int(rows[i][2]) + j*11403 
                
    return(loom)


def spikes_position(trajectory):
    list1 = []
    for j in range(trajectory.number_of_individuals):
        list1 = list1 + [i for i, value in enumerate(trajectory.speed[:,j]) if value > 5]
    return(list1)


def accurate_startles(tr, temp, groupsize, rep):
    list1 = spikes_position(tr)
    loom = loom_frame(temp, groupsize, rep)
    list2 = [value for value in list1 if value < (loom[0] + 1000) and value > (loom[0]) ]
    return(list2) 
    

def first_startle(tr, temp, groupsize, rep):
    a = accurate_startles(tr, temp, groupsize, rep)
    if not a:
        return(accurate_startles(tr, temp, groupsize, rep))
    else:
       return(min(a)) 
    

def latency(tr, temp, groupsize, rep):
    a = first_startle(tr, temp, groupsize, rep)
    if not a:
        return(float('nan'))
    else:
        b = loom_frame(temp, groupsize, rep)
        return(a - b[0])
        
    
    
    
    
    
    
##############################
        
temperature = range(9,30,4)



group = [1,2,4,8,16]



replication = range(3) # number of replicates per treatment


average_latency = np.empty([len(temperature), len(group)])
average_latency.fill(np.nan)


ii = 0 # to keep count of temperature

for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        average_replicate_latency = np.empty([len(replication), 1])
        average_replicate_latency.fill(np.nan)
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
            average_replicate_latency[k] = latency(tr,i,j,k+1)
            
        average_latency[ii,jj] = np.nanmean(average_replicate_latency)
        
        jj= jj + 1
        
    ii = ii + 1