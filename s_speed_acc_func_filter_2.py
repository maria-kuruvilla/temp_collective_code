"""
Created on 10/19/2020

@author: Maria Kuruvilla

Goal - to create own function to calculate speed, accleration and create filter based on distance from origin. (not for individuals)

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
                                 

def track_check_position_masked(tr, temp, group, rep): #replicates start from 1
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for i in range(tr.number_of_individuals):
        
        ax.plot(filter(tr,5)[0], filter(tr,5)[1])
        
    
    ax.set_xlabel('X (BL)')
    ax.set_ylabel('Y (BL)')
    ax.set_title('Trajectories')
    return(ax)

def track_check_masked(tr, temp, group, rep): #replicates start from 1
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    frame_range = range(filter_speed(tr,5).shape[0])
    for i in range(tr.number_of_individuals):
        
        ax.plot(np.asarray(frame_range),filter_speed(tr,5)[frame_range,i])
        
    
    ax.set_xlabel('Frame number')
    ax.set_ylabel('Speed (BL/s)')
    ax.set_title('Temp:' + str(temp) + ' Group:' + str(group) + ' Replicate:' + str(rep))
    return(ax)


def track_check_acc_masked(tr, temp, group, rep): #replicates start from 1
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    frame_range = range(filter_acc(tr,5).shape[0])
    for i in range(tr.number_of_individuals):
        
        ax.plot(np.asarray(frame_range),filter_acc(tr,5)[frame_range,i])
        
    
    ax.set_xlabel('Frame number')
    ax.set_ylabel('Acceleration (BL/s^2)')
    ax.set_title('Temp:' + str(temp) + ' Group:' + str(group) + ' Replicate:' + str(rep))
    return(ax)

