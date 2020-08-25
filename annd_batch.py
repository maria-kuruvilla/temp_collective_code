# -*- coding: utf-8 -*-
"""
Created on Sun May 10 2020

@author: Maria Kuruvilla

Goal - Code to analyse all the tracked videos and calculate annd and save it as pickled file.
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
parser.add_argument('-b', '--integer_b', default=4, type=int)
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
replication = range(args.integer_b) # number of replicates per treatment

#output parent directory
parent_dir = '../../output/temp_collective/'

ii = 0 # to keep count of temperature



def annd(trajectory):
    nnd = np.empty([trajectory.s.shape[0], trajectory.number_of_individuals])
    nnd.fill(np.nan)
    
    nd = np.empty([trajectory.s.shape[0],trajectory.number_of_individuals])
    nd.fill(np.nan)
    
    for i in range(trajectory.number_of_individuals):
        for j in range(trajectory.number_of_individuals):
            if i!=j:
                nd[:,j] = np.sqrt((trajectory.s[:,i,0] - trajectory.s[:,j,0])**2 + (trajectory.s[:,i,1] - trajectory.s[:,i,1])**2)
            
        nnd[:,i] = np.nanmin(nd,1)
        
    annd = np.nanmean(nnd)
        
                
    return(annd)

annd_values = np.empty([len(temperature), len(group)])
annd_values.fill(np.nan)

std_annd_values = np.empty([len(temperature), len(group)])
std_annd_values.fill(np.nan)

for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
        average_replicate_annd = np.empty([len(replication), 1])
        average_replicate_annd.fill(np.nan)
        for k in replication:
                
            input_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
            input_file = input_dir + str(k+1) + '.p'
            tr = pickle.load(open(input_file, 'rb')) # 'rb is for read binary
            average_replicate_annd[k] = annd(tr)
            

        annd_values[ii,jj] = np.nanmean(average_replicate_annd)
        std_annd_values[ii,jj] = np.nanstd(average_replicate_annd)
        jj= jj + 1
        
    ii = ii + 1

out_dir = '../../output/temp_collective/'

# save it as a pickle file
out_fn1 = out_dir + 'annd.p'
pickle.dump(annd_values, open(out_fn1, 'wb')) # 'wb' is for write binary

out_fn2 = out_dir + 'annd_std.p'
pickle.dump(std_annd_values, open(out_fn2, 'wb')) # 'wb' is for write binary