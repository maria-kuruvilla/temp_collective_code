"""
Created on Fri Oct 23  2020

@author: Maria Kuruvilla

Goal - to calculate average and std of many parameters with masked trajectory data
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
    

def local_polarization_func(tr,n=1,x=1): #shape returns tr.s.shape -2 
    if n<x :
        return(np.nan, np.nan)
    if tr.number_of_individuals == 1:
        return(np.nan, np.nan)
    if x >= tr.number_of_individuals:
        return(np.nan, np.nan)
    else:
        indices = ttsocial.neighbour_indices(tr.s[1:-1], n) # indices of the closest n neighbors
        a = ttsocial.restrict(tr.e[1:-1],indices) # normalized velocity (direction) vectors of focal individual and closest n neighbors 
        b = np.multiply(a[:,:,0,:],a[:,:,x,:]) #polarization between focal and xth closest individual
        c=np.sum(b,axis = 2) #dot product of vectors.
        lp = tt.norm(tt.collective.polarization(a))
        #return(c)
        return(np.nanmean(lp[~filter_speed(tr,5).mask]), np.nanmean(c[~filter_speed(tr,5).mask]))
        
def annd(trajectory): #nnd shape returns tr.s.shape -2 
    nnd = np.empty([trajectory.s.shape[0]-2, trajectory.number_of_individuals])
    nnd.fill(np.nan)
    
    nd = np.empty([trajectory.s.shape[0]-2,trajectory.number_of_individuals])
    nd.fill(np.nan)
    
    for i in range(trajectory.number_of_individuals):
        for j in range(trajectory.number_of_individuals):
            if i!=j:
                nd[:,j] = np.sqrt((trajectory.s[1:-1,i,0] - trajectory.s[1:-1,j,0])**2 + (trajectory.s[1:-1,i,1] - trajectory.s[1:-1,i,1])**2)
            
        nnd[:,i] = np.nanmin(nd,1)
        
    annd = np.nanmean(nnd[~filter_speed(tr,5).mask])
        
                
    return(annd)
"""
def spikes(trajectory):
    list1 = []
    for j in range(trajectory.number_of_individuals):
        list1 = list1 + [i for i, value in enumerate(trajectory.speed[:,j]) if value > 10]
    return(len(list1)/trajectory.number_of_individuals)

def startles_total(trajectory):
    list1 = []
    for j in range(trajectory.number_of_individuals):
        list1 = list1 + [i for i, value in enumerate(trajectory.speed[:,j]) if value > 10]
    return(len(list1)/trajectory.number_of_individuals)

def spikes_position(trajectory):
    list1 = []
    for j in range(trajectory.number_of_individuals):
        list1 = list1 + [i for i, value in enumerate(trajectory.speed[:,j]) if value > 10]
    return(list1)
"""
def spikes_position_new(tr): #uses filter_speed
    list1 = []
    for j in range(tr.number_of_individuals):
        list2 = [i for i, value in enumerate(filter_speed(tr,5)[:,j]) if value > 10]
        list2.insert(0,100000000)
        list1 = list1 + [value for i,value in enumerate(list2[1:]) if  (value != (list2[i]+1))]
        
    return(list1)
"""

def get_average_aligment_score(t, number_of_neighbours = 3):
    indices = ttsocial.give_indices(t.s, number_of_neighbours)
    en = ttsocial.restrict(t.e,indices)[...,1:,:]
    alignment = np.nanmean(tt.dot(np.expand_dims(t.e,2), en), axis = -1)
    return np.nanmedian(alignment, axis = -1)
"""

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

average_speed = np.empty([len(temperature), len(group)]) # empty array to calculate average speed per treatment 
average_speed.fill(np.nan)

std_speed = np.empty([len(temperature), len(group)]) # empty array to calculate average speed per treatment 
std_speed.fill(np.nan)

average_acceleration = np.empty([len(temperature), len(group)]) # empty array to calculate average acceleration per treatment
average_acceleration.fill(np.nan)

std_acceleration = np.empty([len(temperature), len(group)]) # empty array to calculate average acceleration per treatment
std_acceleration.fill(np.nan)

annd_values = np.empty([len(temperature), len(group)])
annd_values.fill(np.nan)

std_annd_values = np.empty([len(temperature), len(group)])
std_annd_values.fill(np.nan)
"""
spikes_number = np.empty([len(temperature), len(group)])
spikes_number.fill(np.nan)

std_spikes_number = np.empty([len(temperature), len(group)])
std_spikes_number.fill(np.nan)
"""
polarization = np.empty([len(temperature), len(group)])
polarization.fill(np.nan)

std_polarization = np.empty([len(temperature), len(group)])
std_polarization.fill(np.nan)

#difference_total_accurate = np.empty([len(temperature), len(group)])
#difference_total_accurate.fill(np.nan)

std_difference_total_accurate = np.empty([len(temperature), len(group)])
std_difference_total_accurate.fill(np.nan)

latency_values = np.empty([len(temperature), len(group)])
latency_values.fill(np.nan)

std_latency = np.empty([len(temperature), len(group)])
std_latency.fill(np.nan)

local_polarization = np.empty([len(temperature), len(group)])
local_polarization.fill(np.nan)

std_local_polarization = np.empty([len(temperature), len(group)])
std_local_polarization.fill(np.nan)

ii = 0 # to keep count of temperature

#frames = 5000 #number of frames for which annd is calculated

for i in temperature:
    print(i)
    jj = 0 # to keep count of groups
    for j in group:
        out_dir = parent_dir + '/' + str(i) + '/' + str(j) + '/' 
        
        average_replicate_acceleration = np.empty([len(replication), 1])
        average_replicate_acceleration.fill(np.nan)
            
        average_replicate_speed = np.empty([len(replication), 1])
        average_replicate_speed.fill(np.nan)
        
        average_replicate_annd = np.empty([len(replication), 1])
        average_replicate_annd.fill(np.nan)
        
        #average_replicate_spikes = np.empty([len(replication), 1])
        #average_replicate_spikes.fill(np.nan)
        
        #average_replicate_polarization = np.empty([len(replication), 1])
        #average_replicate_polarization.fill(np.nan)
        
        #difference_total_accurate_replicate = np.empty([len(replication), 1])
        #difference_total_accurate_replicate.fill(np.nan)

        average_replicate_latency = np.empty([len(replication), 1])
        average_replicate_latency.fill(np.nan)
        
        average_replicate_local_polarization = np.empty([len(replication), 1])
        average_replicate_local_polarization.fill(np.nan)
        
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
             
            average_replicate_speed[k] = np.nanmean(filter_speed(tr,5))
            average_replicate_acceleration[k] = np.nanmean(filter_acc(tr,5))
            average_replicate_annd[k] = annd(tr)
            #average_replicate_spikes[k] = spikes(tr)
            #average_replicate_polarization[k] = np.nanmean(tt.norm(tt.collective.polarization(tr.e)))
            #difference_total_accurate_replicate[k] = startles_total(tr) - accurate_startles(tr, i, j, k+1)
            average_replicate_latency[k] = latency(tr,i,j,k+1)
            average_replicate_local_polarization[k] = local_polarization_func(tr,1,1)[1]
            
            
        average_speed[ii,jj] = np.nanmean(average_replicate_speed)
        std_speed[ii,jj] = np.nanstd(average_replicate_speed)
        average_acceleration[ii,jj] = np.nanmean(average_replicate_acceleration)
        std_acceleration[ii,jj] = np.nanstd(average_replicate_acceleration)
        annd_values[ii,jj] = np.nanmean(average_replicate_annd)
        std_annd_values[ii,jj] = np.nanstd(average_replicate_annd)
        #spikes_number[ii,jj] = np.nanmean(average_replicate_spikes)
        #std_spikes_number[ii,jj] = np.nanstd(average_replicate_spikes)
        
        #polarization[ii,jj] = np.nanmean(average_replicate_polarization)
        #std_polarization[ii,jj] = np.nanstd(average_replicate_polarization)
        #difference_total_accurate[ii,jj] = np.nanmean(difference_total_accurate_replicate)
        #std_difference_total_accurate[ii,jj] = np.nanstd(difference_total_accurate_replicate)
        latency_values[ii,jj] = np.nanmean(average_replicate_latency)
        std_latency[ii,jj] = np.nanstd(average_replicate_latency)
        
        local_polarization[ii, jj] = np.nanmean(average_replicate_local_polarization)
        std_local_polarization[ii,jj] = np.nanstd(average_replicate_local_polarization)
        
        jj= jj + 1
        
    ii = ii + 1

out_dir = '../../output/temp_collective/roi/'

# save it as a pickle file
annd_fn1 = out_dir + 'annd.p'
pickle.dump(annd_values, open(annd_fn1, 'wb')) # 'wb' is for write binary

annd_fn2 = out_dir + 'annd_std.p'
pickle.dump(std_annd_values, open(annd_fn2, 'wb')) # 'wb' is for write binary

speed_fn1 = out_dir + 'speed.p'
pickle.dump(average_speed, open(speed_fn1, 'wb')) # 'wb' is for write binary

speed_fn2 = out_dir + 'speed_std.p'
pickle.dump(std_speed, open(speed_fn2, 'wb')) # 'wb' is for write binary

acceleration_fn1 = out_dir + 'acceleration.p'
pickle.dump(average_acceleration, open(acceleration_fn1, 'wb')) # 'wb' is for write binary

acceleration_fn2 = out_dir + 'acceleration_std.p'
pickle.dump(std_acceleration, open(acceleration_fn2, 'wb')) # 'wb' is for write binary
"""
spikes_fn1 = out_dir + 'spikes.p'
pickle.dump(spikes_number, open(spikes_fn1, 'wb')) # 'wb' is for write binary

spikes_fn2 = out_dir + 'spikes_std.p'
pickle.dump(std_spikes_number, open(spikes_fn2, 'wb')) # 'wb' is for write binary

polarization_fn1 = out_dir + 'polarization.p'
pickle.dump(polarization, open(polarization_fn1, 'wb')) # 'wb' is for write binary

polarization_fn2 = out_dir + 'polarization_std.p'
pickle.dump(std_polarization, open(polarization_fn2, 'wb')) # 'wb' is for write binary


accurate_fn1 = out_dir + 'accurate.p'
pickle.dump(difference_total_accurate, open(accurate_fn1, 'wb')) # 'wb' is for write binary

accurate_fn2 = out_dir + 'accurate_std.p'
pickle.dump(std_difference_total_accurate, open(accurate_fn2, 'wb')) # 'wb' is for write binary
"""

latency_fn1 = out_dir + 'latency.p'
pickle.dump(latency_values, open(latency_fn1, 'wb')) # 'wb' is for write binary

latency_fn2 = out_dir + 'latency_std.p'
pickle.dump(std_latency, open(latency_fn2, 'wb')) # 'wb' is for write binary

local_pol_fn1 = out_dir + 'local_pol_m.p'
pickle.dump(local_polarization, open(local_pol_fn1 , 'wb')) # 'wb' is for write binary

local_pol_fn2 = out_dir + 'local_pol_m_std.p'
pickle.dump(std_local_polarization, open(local_pol_fn2, 'wb')) # 'wb' is for write binary

