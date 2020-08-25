# -*- coding: utf-8 -*-
"""
Created on Sun May 10 2020

@author: Maria Kuruvilla

Goal - Code to analyse all the tracked videos and save it as a pickled file.
"""


import sys, os
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
parser.add_argument('-b', '--integer_b', default=3, type=int)
parser.add_argument('-c', '--float_c', default=1.5, type=float)
parser.add_argument('-v', '--verbose', default=True, type=boolean_string)
# Note that you assign a short name and a long name to each argument.
# You can use either when you call the program, but you have to use the
# long name when getting the values back from "args".

# get the arguments
args = parser.parse_args()

#temperatures used in the experiment - used for files naminf=g        
temperature = range(9,30,4)

#group sizes used in the experiment - used for naming files
group = [1,2,4,8,16]

#latest tracked replicate
replication =range(args.integer_b) # number of replicates per treatment

#output parent directory
parent_dir = '../../output/temp_collective/roi/'

ii = 0 # to keep count of temperature


for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir  + str(i) + '/' + str(j) + '/' 
        for k in replication:
            if j==1:
                
                trajectories_file_path = '../../data/temp_collective/roi/'+str(i)+'/' +str(j)+'/GS_'+str(j)+'_T_'+str(i)+'_roi_'+str(k+1)+'/trajectories.npy'
                
            
                    
                
            else:
                
                trajectories_file_path = '../../data/temp_collective/roi/'+str(i)+'/' +str(j)+'/GS_'+str(j)+'_T_'+str(i)+'_roi_'+str(k+1)+'/trajectories_wo_gaps.npy'
            
                
            
            sigma_values = 1.5 #smoothing parameter
            try:
                tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length') # normalizing by body length
            except FileNotFoundError:
                print('File not found')
                print(i,j,k+1)
                continue

            tr.new_time_unit(tr.params['frame_rate'], 'seconds') # changing time unit to seconds
            
            out_fn = out_dir + str(k+1) + '.p'
            pickle.dump(tr, open(out_fn, 'wb')) # 'wb' is for write binary

        jj= jj + 1
        
    ii = ii + 1

