# -*- coding: utf-8 -*-
"""
Created on Mon May 11 2020

@author: Maria Kuruvilla

Goal - Code to filter the data from the tank edges of all the tracked videos and it as pickled file.
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
import pandas as pd

def filter(tr,l):
    left_edge = np.nanmin(tr.s[:,:,0])
    right_edge = np.nanmax(tr.s[:,:,0])
    bottom_edge = np.nanmin(tr.s[:,:,1])
    top_edge = np.nanmax(tr.s[:,:,1])
    column_names = list(range(1,tr.number_of_individuals+1))
    #for i in range(tr.number_of_individuals):
     #   position_x = pd.DataFrame(tr.s[:,i,0], columns = column_names)


    position_x = pd.DataFrame(tr.s[:,:,0], columns = column_names)
    position_y = pd.DataFrame(tr.s[:,:,1], columns = column_names)
    speed = pd.DataFrame(tr.speed[:,:], columns = column_names)
    acceleration = pd.DataFrame(tr.acceleration[:,:], columns = column_names)
    position_x_filtered = position_x.mask((position_x < left_edge + l) | (position_x > right_edge - l))

    position_y_filtered = position_y.mask((position_y < bottom_edge + l) | (position_y > top_edge - l))
    speed_filtered = speed.mask(position_x_filtered.isna() | position_y_filtered.isna())
    acceleration_filtered = acceleration.mask(position_x_filtered.isna() | position_y_filtered.isna())
    x = []
    y = []
    s = []
    a = []
    for i in range(tr.number_of_individuals):
        x.append('x' + str(i+1))
        y.append('y' + str(i+1))
        s.append('speed' + str(i+1))
        a.append('acceleration' + str(i+1))

    position_x_filtered.columns = x
    position_y_filtered.columns = y
    speed_filtered.columns = s
    acceleration_filtered.columns = a

    filtered = pd.concat([position_x_filtered, position_y_filtered, speed_filtered, acceleration_filtered], axis = 1)

    

    return(filtered)

def pandas(tr):

    for i in range(tr.number_of_individuals):
        if i == 0:
            pandas_tr = pd.DataFrame(tr.s[:,i,:], columns = ['x'+str(i) , 'y'+str(i) ])
        else:
            pandas_tr1 = pd.DataFrame(tr.s[:,i,:], columns = ['x'+str(i)  , 'y'+str(i) ])
            pandas_tr = pd.concat([pandas_tr,pandas_tr1], axis = 1)

    return(pandas_tr)

def filtered_track_check(tr, temp, group, rep, l): #replicates start from 1
    frame_range = range(tr.s.shape[0])
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    filtered = filter(tr,l)
    for focal in range(tr.number_of_individuals):
        ax.plot(np.asarray(frame_range),filtered['speed'+str(focal+1)])
        
    
    ax.set_xlabel('Frame number')
    ax.set_ylabel('Speed (BL/s)')
    ax.set_title('Temp:' + str(temp) + ' Group:' + str(group) + ' Replicate:' + str(rep))
    return(ax)



