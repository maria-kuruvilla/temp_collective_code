# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 12:57:41 2020

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

        
temperature = range(9,30,4)



group = [1,2,4,8,16]


polarization = np.array([[1.        , 0.74083665, 0.62935441, 0.4908591 , 0.41110361],
       [1.        , 0.75069572, 0.59222193, 0.46613094, 0.48811376],
       [1.        , 0.74265335, 0.69532405, 0.54092712, 0.34672795],
       [1.        , 0.75540608, 0.59358043, 0.43391956, 0.35321903],
       [1.        , 0.67285949, 0.6224674 , 0.50056763, 0.32031768],
       [1.        , 0.80708511, 0.66953274, 0.49475468, 0.31663836]])

std_polarization = np.array([[0.        , 0.06128706, 0.00478066, 0.02072101, 0.08162972],
       [0.        , 0.08303297, 0.03176809, 0.03654017, 0.06559605],
       [0.        , 0.03271406, 0.06417445, 0.08609299, 0.04542871],
       [0.        , 0.04117595, 0.03278339, 0.02294503, 0.00973238],
       [0.        , 0.0194057 , 0.11449309, 0.03979405, 0.01622424],
       [0.        , 0.03427955, 0.03426918, 0.0509842 , 0.00184274]])

for i in range(6):
    plt.plot(group, polarization[i,:], label = str(temperature[i]), linewidth = 0.5)
    plt.fill_between(group, polarization[i,:] - std_polarization[i,:], polarization[i,:] + std_polarization[i,:], alpha = 0.3)
    
plt.xlabel('Group Size')
plt.ylabel('Polarization')
plt.legend()