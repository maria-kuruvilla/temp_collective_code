"""
Goal - to enter the temp, group size and replicate and func should return plot with speed as a function of time
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

def spikes_position(trajectory):
    list1 = []
    for j in range(trajectory.number_of_individuals):
        list2 = [i for i, value in enumerate(trajectory.speed[:,j]) if value > 10]
        list2.insert(0,100000000)
        list1 = list1 + [value for i,value in enumerate(list2[1:]) if  (value != (list2[i]+1))]
        
    return(list1)

def track_check(tr, temp, group, rep): #replicates start from 1
    frame_range = range(tr.s.shape[0])
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for focal in range(tr.number_of_individuals):
        ax.plot(np.asarray(frame_range),tr.speed[frame_range, focal])
        
    
    ax.set_xlabel('Frame number')
    ax.set_ylabel('Speed (BL/s)')
    ax.set_title('Temp:' + str(temp) + ' Group:' + str(group) + ' Replicate:' + str(rep))
    return(ax)

def track_check_acc(tr, temp, group, rep): #replicates start from 1
    frame_range = range(tr.s.shape[0])
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for focal in range(tr.number_of_individuals):
        ax.plot(np.asarray(frame_range),tr.acceleration[frame_range, focal])
        
    
    ax.set_xlabel('Frame number')
    ax.set_ylabel('Acceleration (BL/s^2)')
    ax.set_title('Temp:' + str(temp) + ' Group:' + str(group) + ' Replicate:' + str(rep))
    return(ax)

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
input_file = input_dir + str(args.integer_b3) + '.p'
sigma_values = 1.5 #smoothing parameter
try:
    tr = pickle.load(open(input_file, 'rb')) # 'rb is for read binary
except FileNotFoundError:
    print(args.integer_b1,args.integer_b2,args.integer_b3)
    print('File not found')
    pass
track_check(tr, args.integer_b1, args.integer_b2, args.integer_b3)
track_check_acc(tr,args.integer_b1,args.integer_b2,args.integer_b3)
print(spikes_position(tr))
plt.show()