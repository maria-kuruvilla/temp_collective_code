"""
Goal - to produce a csv file with temp, gs rep and 90th percentile speed to do stats
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

def spikes_position_new(tr): #uses filter_speed
    list1 = []
    for j in range(tr.number_of_individuals):
        list2 = [i for i, value in enumerate(filter_speed(tr,5)[:,j]) if value > 10]
        list2.insert(0,100000000)
        list1 = list1 + [value for i,value in enumerate(list2[1:]) if  (value != (list2[i]+1))]
        
    return(list1)

rows = []
with open('../../data/temp_collective/looms_roi.csv', 'r') as csvfile:
    looms = csv.reader(csvfile)
    for row in looms:
        rows.append(row)
        
        


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

temperature = range(9,30,4)



group = [1,2,4,8,16,32]



replication = range(10) # number of replicates per treatment

#output parent directory
parent_dir = '../../output/temp_collective/roi'

with open('../../data/temp_collective/roi/stats_speed_acc_latency.csv', mode='w') as stats_speed:
    writer = csv.writer(stats_speed, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    writer.writerow(['Temperature', 'Groupsize', 'Replicate', '90_speed','90_acc','latency'])
    for i in temperature:
        print(i)
        jj = 0 # to keep count of groups
        for j in group:
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
                perc_speed = np.percentile(filter_speed(tr,5).compressed(),90)
                perc_acc = np.percentile(filter_acc(tr,5).compressed(),90)
                lat = latency(tr,i,j,k+1)
                writer.writerow([i, j, k+1,perc_speed,perc_acc,lat])

