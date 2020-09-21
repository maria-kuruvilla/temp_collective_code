"""
Goal - to enter the temp, group size and replicate and func should return a histogram of speed
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


def speed_histogram(x): 
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    #speed_pdf, speed = np.histogram(x,bins=10,density=True)
    #plt.plot(speed,np.log(speed_pdf))
    ax.hist(x, density=True, bins=30, log = True)
    #ax.set_xscale('log')
    #ax.set_xlim(left = 5)
    #ax.set_ylim([0,0.0002])
    ax.set_xlabel('Speed (BL/s)')
    ax.set_ylabel('Probability')
    #plt.show()
    #plt.xticks(ticks = [10,20,100,300], labels = [10,20,100,300])
    
    out_dir = parent_dir = '../../output/temp_collective/roi_figures/speed_pdf_log_lin_no_smooth.png'
    fig.savefig(out_dir, dpi = 300)
    return(ax)
    

def acc_histogram(y): #replicates start from 1
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    ax.hist(y, density=True, bins=100, log = True)
    #ax.set_xscale('log')
    #plt.xticks(ticks = [10,100,500,1000,2000,5000], labels = [10,100,500,1000,2000,5000])
    #ax.set_xlim(left = 5)
    #ax.set_ylim([0,0.0002])
    ax.set_xlabel('Acceleration')
    ax.set_ylabel('Probability')
    #ax.set_title('Temp:' + str(temp) + ' Group:' + str(group) + ' Replicate:' + str(rep))
    out_dir = parent_dir = '../../output/temp_collective/roi_figures/acc_pdf_log_lin_no_smooth.png'
    fig.savefig(out_dir, dpi = 300)
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
parser.add_argument('-b', '--integer_b', default=10, type=int)
parser.add_argument('-c', '--float_c', default=1.5, type=float)
parser.add_argument('-v', '--verbose', default=True, type=boolean_string)
# Note that you assign a short name and a long name to each argument.
# You can use either when you call the program, but you have to use the
# long name when getting the values back from "args".

# get the arguments
args = parser.parse_args()

temperature = range(9,30,4)



group = [1,2,4,8,16]



replication = range(args.integer_b) # number of replicates per treatment

#output parent directory
parent_dir = '../../output/temp_collective/roi'

ii = 0 # to keep count of temperature
"""
#frames = 5000 #number of frames for which annd is calculated
x = [] #append speed values to this 
y=[]
for i in temperature:
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
        for k in replication:
            
            input_file = out_dir + str(k+1) + '_nosmooth.p'
            
            try:
                tr = pickle.load(open(input_file, 'rb')) # 'rb is for read binary
            except FileNotFoundError:
                print(i,j,k)
                print('File not found')
                continue
            #for m in range(tr.speed.shape[1]):
             #   x= np.r_[x,tr.speed[:,m]]
            for m in range(tr.acceleration.shape[1]):
                y= np.r_[y,tr.acceleration[:,m]]
                x= np.r_[x,tr.speed[:,m]]
"""
acc_histogram(y)
speed_histogram(x)
#acc_histogram(tr, args.integer_b1, args.integer_b2, args.integer_b3)
plt.show()