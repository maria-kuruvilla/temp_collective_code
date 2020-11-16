"""
Created on Fri Nov 6 2020

@author: Maria Kuruvilla

Goal - number of startles total per unit unmasked time (per hour)
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
parser.add_argument('-b', '--integer_b', default=10, type=int)
parser.add_argument('-c', '--float_c', default=1.5, type=float)
parser.add_argument('-v', '--verbose', default=True, type=boolean_string)
# Note that you assign a short name and a long name to each argument.
# You can use either when you call the program, but you have to use the
# long name when getting the values back from "args".

# get the arguments
args = parser.parse_args()

#functions
    
def position(tr): #shape returns tr.s.shape
    return(tr.s)


def speed(tr): #speed(tr).shape returns tr.speed.shape - 2
    v = (position(tr)[2:] - position(tr)[:-2]) / 2
    b = np.linalg.norm(v, axis=-1)
    return(b*60)

def acceleration(tr): #shape returns tr.acceleration.shape - 2
    a = position(tr)[2:] - 2 * position(tr)[1:-1] + position(tr)[:-2]
    aa = np.linalg.norm(a, axis=-1)  
    return(aa*3600)
        
def e(tr): #e.shape returns tr.speed.shape - 2
    vel = (position(tr)[2:] - position(tr)[:-2]) / 2
    n = np.linalg.norm(v,axis = 2)  
    return(vel/n[...,np.newaxis])

    

def filter(tr, roi = 5): #ind (for individual) starts from 0, roi - edge of region of interest
#shape returns tr.s.shape -2 
    position_mask0 = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), position(tr)[1:-1,:,0],copy=False)
    position_mask1 = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), position(tr)[1:-1,:,1],copy=False)
    return(np.ma.dstack((position_mask0,position_mask1)))  
    
def filter_speed(tr, roi = 5): #shape returns tr.s.shape -2 
    speed_mask = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), speed(tr),copy=False)
    
    return(speed_mask)         



def filter_acc(tr, roi = 5): #shape returns tr.s.shape -2 
    acc_mask = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), acceleration(tr),copy=False)
    
    return(acc_mask)#[~acc_mask.mask].data) 

    
def filter_e(tr, roi =5):
    x = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), e(tr)[:,:,0],copy=False)
    y = np.ma.masked_where((tr.distance_to_origin[1:-1] > roi)|(tr.distance_to_origin[0:-2] > roi)|(tr.distance_to_origin[2:] > roi), e(tr)[:,:,1],copy=False)
    return(np.ma.dstack((x,y)))  
    
def spikes_position_new(tr): #uses filter_speed
    list1 = []
    for j in range(tr.number_of_individuals):
        list2 = [i for i, value in enumerate(filter_speed(tr,5)[:,j]) if value > 10]
        list2.insert(0,100000000)
        list1 = list1 + [value for i,value in enumerate(list2[1:]) if  (value != (list2[i]+1))]
        
    return(list1)

def startles_list_new(tr): #uses filtered data
    return(len(spikes_position_new(tr)))

rows = []
with open('../../data/temp_collective/looms_roi.csv', 'r') as csvfile:
    looms = csv.reader(csvfile)
    for row in looms:
        rows.append(row)
        
        

        
temperature = range(9,30,4)



group = [1,2,4,8,16,32]



replication = range(args.integer_b) # number of replicates per treatment

#output parent directory
parent_dir = '../../output/temp_collective/roi'

def loom_frame(temp, groupsize, rep):
    if temp == 29:
        cam = 'Cam 7'
    elif temp == 25:
        cam = 'Cam 8'
    elif temp == 17:
        cam = 'Cam 9'
    elif temp == 13:
        cam = 'Cam 10'
    elif temp == 21:
        cam = 'Cam 11'
    elif temp == 9:
        cam = 'Cam 12'
    g = str(groupsize)
    r = str(rep)
    loom = np.zeros([5,1])        
    for i in range(len(rows)):
        if rows[i][1]==cam and rows[i][3]==g and rows[i][4]==r:
            for j in range(5):
                loom[j] = int(rows[i][2]) + j*11403 
    
    return(loom)

def accurate_startles(tr, temp, groupsize, rep): #uses filtered speed
    list1 = spikes_position_new(tr)
    loom = loom_frame(temp, groupsize, rep)
    list2 = [i for i, value in enumerate(list1[:]) if value < (loom[0] + 700) and value > (loom[0]+500) ]
    list2 = list2 + [i for i, value in enumerate(list1[:]) if value < (loom[1] + 700) and value > (loom[1]+500) ]
    list2 = list2 + [i for i, value in enumerate(list1[:]) if value < (loom[2] + 700) and value > (loom[2]+500) ]
    list2 = list2 + [i for i, value in enumerate(list1[:]) if value < (loom[3] + 700) and value > (loom[3]+500) ]
    list2 = list2 + [i for i, value in enumerate(list1[:]) if value < (loom[4] + 700) and value > (loom[4]+500) ]
    
    return(len(list2)/tr.number_of_individuals)

def accurate_startles_frame(tr, temp, groupsize, rep,i): #i starts from 0 #uses filtered data
    list1 = spikes_position_new(tr)
    loom = loom_frame(temp, groupsize, rep)
    list2 = [value for value in list1 if (value < (loom[i] + 700) and value > (loom[i]+500)) ]
    return(list2) 
    

def first_startle(tr, temp, groupsize, rep,i): #uses filtered data
    a = accurate_startles_frame(tr, temp, groupsize, rep,i) # i starts from 0
    if not a:
        return(np.nan)
    else:
       return(min(a)) 
    

def latency(tr, temp, groupsize, rep): #uses filtred data
    a = np.empty([5,1])
    a.fill(np.nan)
    b = loom_frame(temp, groupsize, rep)

    for i in range(5):
        a[i] = first_startle(tr, temp, groupsize, rep,i) - b[i]
    
    return(np.nanmean(a))

average_startles = np.empty([len(temperature), len(group)]) # empty array to calculate average no. of startles per treatment
average_startles.fill(np.nan)

std_startles = np.empty([len(temperature), len(group)]) # empty array to calculate std of no. of startles per treatment 
std_startles.fill(np.nan)

ii = 0 # to keep count of temperature


for i in temperature:
    print(i)
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
        
        replicate_startles = np.empty([len(replication), 1])
        replicate_startles.fill(np.nan)
        for k in replication:

            if j == 1:
                trajectories_file_path = '../../data/temp_collective/roi/'+str(i)+'/' +str(j)+'/GS_'+str(j)+'_T_'+str(i)+'_roi_'+str(k+1)+'/trajectories.npy'
            else:
                trajectories_file_path = '../../data/temp_collective/roi/'+str(i)+'/' +str(j)+'/GS_'+str(j)+'_T_'+str(i)+'_roi_'+str(k+1)+'/trajectories_wo_gaps.npy'
            try:
                tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, 		           center=True).normalise_by('body_length')
                tr.new_time_unit(tr.params['frame_rate'], 'seconds')
            except FileNotFoundError:
                print(i,j,k)
                print('File not found')
                continue
            
            if  filter_speed(tr).compressed().shape[0] == 0:
                replicate_startles[k] = np.nan
            else:

                replicate_startles[k] = startles_list_new(tr)*60*60*60/filter_speed(tr).compressed().shape[0]

        average_startles[ii, jj] = np.nanmean(replicate_startles)
        std_startles[ii,jj] = np.nanstd(replicate_startles)
        
        jj= jj + 1
        
    ii = ii + 1

out_dir = '../../output/temp_collective/roi/'

# save it as a pickle file
startles_fn1 = out_dir + 'unmasked_startles.p'
pickle.dump(average_startles, open(startles_fn1, 'wb')) # 'wb' is for write binary

startles_fn2 = out_dir + 'unmasked_startles_std.p'
pickle.dump(std_startles, open(startles_fn2, 'wb')) # 'wb' is for write binary
