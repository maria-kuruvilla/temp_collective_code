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
    


    position_x = pd.DataFrame(tr.s[:,:,0])
    position_y = pd.DataFrame(tr.s[:,:,1])

    position_x_filtered = position_x.mask((position_x < left_edge + l) | (position_x > right_edge - l))
    position_y_filtered = position_y.mask((position_y < bottom_edge + l) | (position_y > top_edge - l))


    

    return(tr)

    def pandas(tr):

        for i in range(tr.number_of_individuals):
            if i == 0:
                pandas_tr = pd.DataFrame(tr.s[:,i,:])
            else:
                pandas_tr1 = pd.DataFrame(tr.s[:,i,:])
                pandas_tr = pd.concat([pandas_tr,pandas_tr1], axis = 1)

    return(pandas_tr)



