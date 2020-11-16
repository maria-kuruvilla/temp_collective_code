
"""
Goal - to check if curvature can be used to characterize startle
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
                                 
def filter_curve(tr, roi = 5): 
    curve_mask = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), tr.curvature[1:-1],copy=False)
    
    return(curve_mask)  


def curve(tr,frame0,frame1):
    fig, ax = plt.subplots(2,figsize=(15,15))

    frame_range = range(frame0, frame1,1) 

    for focal in range(tr.number_of_individuals):
        
        ax[0].plot(np.asarray(frame_range),filter_speed(tr,5)[frame_range, focal])
    
        ax[1].plot(np.asarray(frame_range),filter_curve(tr,5)[frame_range, focal])
        plt.subplots_adjust(hspace = 1.)

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
parser.add_argument('-b1', '--integer_b1', default=29, type=int)
parser.add_argument('-b2', '--integer_b2', default=16, type=int)
parser.add_argument('-b3', '--integer_b3', default=3, type=int)
parser.add_argument('-c', '--float_c', default=1.5, type=float)
parser.add_argument('-v', '--verbose', default=True, type=boolean_string)
# Note that you assign a short name and a long name to each argument.
# You can use either when you call the program, but you have to use the
# long name when getting the values back from "args".

# get the arguments
args = parser.parse_args()

parent_dir = '../../output/temp_collective/roi'

input_dir = parent_dir + '/' + str(args.integer_b1) + '/' + str(args.integer_b2) + '/' 
input_file = input_dir + str(args.integer_b3) + '_nosmooth.p'
#sigma_values = 1.5 #smoothing parameter
if args.integer_b2 == 1:
    trajectories_file_path = '../../data/temp_collective/roi/'+str(args.integer_b1)+'/' +str(args.integer_b2)+'/GS_'+str(args.integer_b2)+'_T_'+str(args.integer_b1)+'_roi_'+str(args.integer_b3)+'/trajectories.npy'
else:
    trajectories_file_path = '../../data/temp_collective/roi/'+str(args.integer_b1)+'/' +str(args.integer_b2)+'/GS_'+str(args.integer_b2)+'_T_'+str(args.integer_b1)+'_roi_'+str(args.integer_b3)+'/trajectories_wo_gaps.npy'
try:
    tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, 		                    center=True).normalise_by('body_length')
    tr.new_time_unit(tr.params['frame_rate'], 'seconds')
except FileNotFoundError:
    print(args.integer_b1,args.integer_b2,args.integer_b3)
    print('File not found')
    pass
curve(tr,0,tr.s.shape[0]-2)
plt.show()
