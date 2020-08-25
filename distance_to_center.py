# -*- coding: utf-8 -*-
"""
Created on Aug 20 2020
Goal - To calculate the average distance to center for each treatment

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
    

temperature = range(9,30,4)



group = [1,2,4,8,16]



replication = range(10) # number of replicates per treatment


distance_to_center = np.empty([len(temperature), len(group)])
distance_to_center.fill(np.nan)

std_distance_to_center = np.empty([len(temperature), len(group)])
std_distance_to_center.fill(np.nan)




#output parent directory
parent_dir = '../../output/temp_collective/roi'

ii = 0 # to keep count of temperature

#frames = 5000 #number of frames for which annd is calculated

for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
        
        average_replicate_dtc = np.empty([len(replication), 1])
        average_replicate_dtc.fill(np.nan)

        for k in replication:
            
            input_file = out_dir + str(k+1) + '.p'
            
            try:
                tr = pickle.load(open(input_file, 'rb')) # 'rb is for read binary
            except FileNotFoundError:
                print(i,j,k+1)
                print('File not found')
                continue
             
            average_replicate_dtc[k] = np.nanmean(tr.distance_to_center)
            
            
            
        
        distance_to_center[ii, jj] = np.nanmean(average_replicate_dtc)
        std_distance_to_center[ii,jj] = np.nanstd(average_replicate_dtc)

        
        
        jj= jj + 1
        
    ii = ii + 1

out_dir = '../../output/temp_collective/roi/'

# save it as a pickle file
dtc_fn1 = out_dir + 'dtc.p'
pickle.dump(distance_to_center, open(dtc_fn1, 'wb')) # 'wb' is for write binary

dtc_fn2 = out_dir + 'dtc_std.p'
pickle.dump(std_distance_to_center, open(dtc_fn2, 'wb')) # 'wb' is for write binary