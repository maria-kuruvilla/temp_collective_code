# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 2020

@author: Maria Kuruvilla

Goal - Code to visualise heatmap of position of fish
"""




import os
from pprint import pprint
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy import stats
# trajectorytools needs to be installed. To installed follow the instructions 
# at http://www.github.com/fjhheras/trajectorytools
import trajectorytools as tt
import trajectorytools.plot as ttplot
import trajectorytools.socialcontext as ttsocial
import argparse
import pickle

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



hist2d = []
vmin = 0
vmax = 0
"""
for focal in range(tr.number_of_individuals):
    hist2d.append(np.histogram2d(tr.s[:, focal, 0][~np.isnan(tr.s[:,focal,0])], tr.s[:, focal, 1][~np.isnan(tr.s[:,focal,1])], 25)[0])
    vmax = max(vmax, hist2d[focal].max()) 
"""
hist2d.append(np.histogram2d(tr.s[:, :, 0][~np.isnan(tr.s[:,:,0])], tr.s[:, :, 1][~np.isnan(tr.s[:,:,1])], 12)[0])
vmax = np.max(hist2d) 

# Plot distributions of positions in the arena
figv, ax_hist = plt.subplots(1,figsize=(15,7), sharey=True, sharex=True)
min_x, max_x = np.nanmin(tr.s[..., 0]), np.nanmax(tr.s[..., 0])
min_y, max_y = np.nanmin(tr.s[..., 1]), np.nanmax(tr.s[..., 1])
"""
for focal in range(tr.number_of_individuals):
    ax = ax_hist[focal]
    ax.imshow(np.rot90(hist2d[focal]), interpolation='none',
              extent=[min_x, max_x, min_y, max_y],
              vmin=vmin, vmax=vmax)
    ax.set_title('fish {}'.format(focal), fontsize=14)
    if ax.is_last_row():
            ax.set_xlabel('X position (BL)')
    if ax.is_first_col():
        ax.set_ylabel('Y position (BL)')
    ax.tick_params(labelsize=14)
    ax.set_title("Fish {}".format(focal+1),fontsize=14)
    ax.set_aspect('equal','box')
"""
ax_hist.imshow(np.rot90(hist2d)[:,0,:], interpolation = 'none',extent=[min_x, max_x, min_y, max_y],vmin=vmin, vmax=vmax)
plt.show()