"""
Goal - To calculate x percentile speed for each group size and temperature.

Created - 10/20/20
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

def position(tr):
    return(tr.s)


def speed(tr):
    v = (position(tr)[2:] - position(tr)[:-2]) / 2
    b = np.linalg.norm(v, axis=-1)
    return(b*60)

def acceleration(tr):
    a = position(tr)[2:] - 2 * position(tr)[1:-1] + position(tr)[:-2]
    aa = np.linalg.norm(a, axis=-1)  
    return(aa*3600)
        

def filter(tr, roi = 5): #ind (for individual) starts from 0, roi - edge of region of interest
    position_mask0 = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), position(tr)[1:-1,:,0],copy=False)
    position_mask1 = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), position(tr)[1:-1,:,1],copy=False)
    return(position_mask0,position_mask1)  
    
def filter_speed(tr, roi = 5): 
    speed_mask = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), speed(tr),copy=False)
    
    return(speed_mask)         



def filter_acc(tr, roi = 5): 
    acc_mask = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), acceleration(tr),copy=False)
    
    return(acc_mask)#[~acc_mask.mask].data)       

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
parser.add_argument('-b', '--integer_b', default=90, type=int)
parser.add_argument('-c', '--float_c', default=1.5, type=float)
parser.add_argument('-d', '--integer_d', default=1, type=int)
parser.add_argument('-v', '--verbose', default=True, type=boolean_string)
# Note that you assign a short name and a long name to each argument.
# You can use either when you call the program, but you have to use the
# long name when getting the values back from "args".

# get the arguments
args = parser.parse_args()

#functions
    

temperature = range(9,30,4)



group = [1,2,4,8,16]



replication = range(10) # number of replicates per treatment


percentile_speed = np.empty([len(temperature), len(group)])
percentile_speed.fill(np.nan)

std_percentile_speed = np.empty([len(temperature), len(group)])
std_percentile_speed.fill(np.nan)



#output parent directory
parent_dir = '../../data/temp_collective/roi'

ii = 0 # to keep count of temperature

#frames = 5000 #number of frames for which annd is calculated

for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
        
        average_replicate_percentile_speed = np.empty([len(replication), 1])
        average_replicate_percentile_speed.fill(np.nan)

        for k in replication:
            
            if j == 1:
                input_file = out_dir + '/GS_'+str(j)+'_T_'+str(i)+'_roi_'+str(k+1)+'/trajectories.npy'
            else:   
                input_file = out_dir + '/GS_'+str(j)+'_T_'+str(i)+'_roi_'+str(k+1)+'/trajectories_wo_gaps.npy'
                
            
            
            
            try:
                tr = tt.Trajectories.from_idtrackerai(input_file, 		           center=True).normalise_by('body_length')
                tr.new_time_unit(tr.params['frame_rate'], 'seconds')		
            
            except FileNotFoundError:
                print(i,j,k)
                print('File not found')
                continue
             
            average_replicate_percentile_speed[k] = np.percentile(filter_speed(tr,5).compressed(),args.integer_b)
            
            
            
        
        percentile_speed[ii, jj] = np.nanmean(average_replicate_percentile_speed)
        std_percentile_speed[ii,jj] = np.nanstd(average_replicate_percentile_speed)

        
        
        jj= jj + 1
        
    ii = ii + 1

out_dir = '../../output/temp_collective/roi/'

# save it as a pickle file
p_c_fn1 = out_dir + 'percentile_speed' + str(args.integer_b)+ '.p'
pickle.dump(percentile_speed, open(p_c_fn1, 'wb')) # 'wb' is for write binary

p_c_fn2 = out_dir + 'percentile_speed' + str(args.integer_b) + '_std.p'
pickle.dump(std_percentile_speed, open(p_c_fn2, 'wb')) # 'wb' is for write binary
 
