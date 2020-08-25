# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 2020

@author: Maria Kuruvilla

Goal - Code to visualise trajectoies in the tank to know the dimensions of the tank
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

parent_dir = '../../output/temp_collective/roi'

def track_tank(i,j,k):


	input_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
	input_file = input_dir + str(k) + '.p'
	sigma_values = 1.5 #smoothing parameter
	try:
	    tr = pickle.load(open(input_file, 'rb')) # 'rb is for read binary
	except FileNotFoundError:
	    print(i,j,k)
	    print('File not found')
	    pass

	fig, ax_trajectories = plt.subplots(figsize=(5,5))
	#time_range= (0, 60) # SET HERE THE RANGE IN SECONDS FOR WHICH YOU WANT TO PLOT THE POSITIONS
	frame_range = range(tr.s.shape[0]) 

	for i in range(tr.number_of_individuals):
		ax_trajectories.plot(tr.s[frame_range,i,0], tr.s[frame_range,i,1])
		ax_trajectories.set_aspect('equal','box')
		ax_trajectories.set_title('Trajectories',fontsize=24)
		ax_trajectories.set_xlabel('X (BL)',fontsize=24)
		ax_trajectories.set_ylabel('Y (BL)',fontsize=24)
	plt.show()
"""
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
#print(spikes_position(tr))
plt.show()
"""