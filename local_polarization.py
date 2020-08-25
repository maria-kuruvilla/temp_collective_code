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
import pickle
import argparse

#argparse
def boolean_string(s):
    # this function helps with getting Boolean input
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

# create the parser object
parser = argparse.ArgumentParser()

# NOTE: argparse will throw an error if:
#     - a flag is given with no value
#     - the value does not match the type
# and if a flag is not given it will be filled with the default.
parser.add_argument('-a', '--a_string', default='hi', type=str)
parser.add_argument('-b', '--integer_b', default=11, type=int)
parser.add_argument('-c', '--float_c', default=1.5, type=float)
parser.add_argument('-d', '--integer_d', default=1, type=int)
parser.add_argument('-v', '--verbose', default=True, type=boolean_string)
# Note that you assign a short name and a long name to each argument.
# You can use either when you call the program, but you have to use the
# long name when getting the values back from "args".

# get the arguments
args = parser.parse_args()

#functions
    

def get_average_local_polarization(t, number_of_neighbours = 1):
    
    if t.number_of_individuals == 1:
        return('nan')
    else:
        indices = ttsocial.give_indices(t.s, number_of_neighbours)
        en = ttsocial.restrict(t.e,indices)
        local_polarization = tt.norm(tt.collective.polarization(en))
        #return np.nanmean(local_polarization, axis = -1) 
        return np.nanmean(local_polarization)

def trajectory(i, j , k): #takes replicates starting from 1
    if j == 1:
        trajectories_file_path = '../../data/temp_collective/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k)+'/trajectories.npy'
    else:
        trajectories_file_path = '../../data/temp_collective/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k)+'/trajectories_wo_gaps.npy'
    sigma_values = 1.5 #smoothing parameter
    tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length') # normalizing by body length
    tr.new_time_unit(tr.params['frame_rate'], 'seconds') # changing time unit to seconds
    return(tr) 

def local_polarization(tr,n=1,x=1):
    if n<x :
        return(np.nan, np.nan)
    if tr.number_of_individuals == 1:
        return(np.nan, np.nan)
    if x >= tr.number_of_individuals:
        return(np.nan, np.nan)
    else:
        indices = ttsocial.give_indices(tr.s, n) # indices of the closest n neighbors
        a = ttsocial.restrict(tr.e,indices) # normalized velocity (direction) vectors of focal individual and closest n neighbors 
        b = np.multiply(a[:,:,0,:],a[:,:,x,:]) #polarization between focal and xth closest individual
        c=np.sum(b,axis = 2) #dot product of vectors.
        lp = tt.norm(tt.collective.polarization(a))
        #return(c)
        return(np.nanmean(lp), np.nanmean(c))
        
def global_polarization(tr):
    polarization_order_parameter = tt.norm(tt.collective.polarization(tr.e))
    return(polarization_order_parameter)

temperature = range(9,30,4)



group = [1,2,4,8,16]



replication = range(args.integer_b) # number of replicates per treatment


local_pol = np.empty([len(temperature), len(group)])
local_pol.fill(np.nan)

std_local_pol = np.empty([len(temperature), len(group)])
std_local_pol.fill(np.nan)

local_pol_m = np.empty([len(temperature), len(group)])
local_pol_m.fill(np.nan)

std_local_pol_m = np.empty([len(temperature), len(group)])
std_local_pol_m.fill(np.nan)



#output parent directory
parent_dir = '../../output/temp_collective/roi'

ii = 0 # to keep count of temperature

#frames = 5000 #number of frames for which annd is calculated

for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
        
        average_replicate_local_pol = np.empty([len(replication), 1])
        average_replicate_local_pol.fill(np.nan)

        average_replicate_local_pol_m = np.empty([len(replication), 1])
        average_replicate_local_pol_m.fill(np.nan)


        for k in replication:
            
            input_file = out_dir + str(k+1) + '.p'
            
            try:
                tr = pickle.load(open(input_file, 'rb')) # 'rb is for read binary
            except FileNotFoundError:
                print(i,j,k)
                print('File not found')
                continue
             
            average_replicate_local_pol[k] = local_polarization(tr,1,1)[0]
            average_replicate_local_pol_m[k] = local_polarization(tr,1,1)[1]
            
            
        
        local_pol[ii, jj] = np.nanmean(average_replicate_local_pol)
        std_local_pol[ii,jj] = np.nanstd(average_replicate_local_pol)

        local_pol_m[ii, jj] = np.nanmean(average_replicate_local_pol_m)
        std_local_pol_m[ii,jj] = np.nanstd(average_replicate_local_pol_m)
        
        jj= jj + 1
        
    ii = ii + 1

out_dir = '../../output/temp_collective/roi/'

# save it as a pickle file
local_pol_fn1 = out_dir + 'local_pol.p'
pickle.dump(local_pol, open(local_pol_fn1, 'wb')) # 'wb' is for write binary

local_pol_fn2 = out_dir + 'local_pol_std.p'
pickle.dump(std_local_pol, open(local_pol_fn2, 'wb')) # 'wb' is for write binary

local_pol_m_fn1 = out_dir + 'local_pol_m' + str(x)+ '.p'
pickle.dump(local_pol_m, open(local_pol_m_fn1, 'wb')) # 'wb' is for write binary

local_pol_m_fn2 = out_dir + 'local_pol_m_std.p'
pickle.dump(std_local_pol_m, open(local_pol_m_fn2, 'wb')) # 'wb' is for write binary