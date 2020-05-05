# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 17:19:16 2020

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

#functions
    

def get_average_local_polarization(t, number_of_neighbours = 1):
    t.new_length_unit(t.params['radius_px'], 'radius')
    if t.number_of_individuals == 1:
        return('nan')
    else:
        indices = ttsocial.give_indices(t.s, number_of_neighbours)
        en = ttsocial.restrict(t.e,indices)[...,1:,:]
        local_polarization = tt.norm(tt.collective.polarization(en))
        return np.nanmean(local_polarization, axis = -1) 
    
def trajectory(i, j , k): #takes replicates starting from 1
    if j == 1:
        trajectories_file_path = '../../data/temp_collective/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k)+'/trajectories.npy'
    else:
        trajectories_file_path = '../../data/temp_collective/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k)+'/trajectories_wo_gaps.npy'
    sigma_values = 1.5 #smoothing parameter
    tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length') # normalizing by body length
    tr.new_time_unit(tr.params['frame_rate'], 'seconds') # changing time unit to seconds
    return(tr) 