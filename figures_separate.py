"""
Code to plot all parameters as a function of group size and temperature - separate figures for each group size 
"""

# imports
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm

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
parser.add_argument('-a', '--a_string', default='annd', type=str)
#parser.add_argument('-s', '--a_string', default='annd_std', type=str)
parser.add_argument('-b', '--integer_b', default=3, type=int)
parser.add_argument('-c', '--float_c', default=1.5, type=float)
parser.add_argument('-v', '--verbose', default=True, type=boolean_string)
# Note that you assign a short name and a long name to each argument.
# You can use either when you call the program, but you have to use the
# long name when getting the values back from "args".

# get the arguments
args = parser.parse_args()
xx=0
h = 0.3
if args.a_string=='annd':
	y_label = 'ANND (Body Length)'
	xx = 1
if args.a_string=='speed':
	y_label = 'Speed (Body Length/s)'
if args.a_string=='acceleration':
	y_label = 'Acceleration (Body Length/s'+r'$^2$)'
if args.a_string=='polarization':
	y_label = 'Polarization'
	xx=1
if args.a_string=='spikes':
	y_label = 'Number of \n startles'
	h = 0.4
if args.a_string=='accurate':
	y_label = 'Number of \n accurate startles'
	h = 0.4 
if args.a_string=='latency':
	y_label = 'Latency (frames)'
if args.a_string=='local_pol':
	y_label = 'Local polarization'
	xx = 1
if args.a_string=='local_pol_m':
	y_label = 'Local polarization'
	xx = 1
if args.a_string=='dtc':
	y_label = 'Distance to center \n (pixels)'

if args.a_string=='percentile_speed99':
	y_label = '99th percentile of speed \n (Body Length/s)'

if args.a_string=='percentile_speed90':
	y_label = '90th percentile of speed \n (Body Length/s)'
if args.a_string=='percentile_speed80':
	y_label = '80th percentile of speed \n (Body Length/s)'
if args.a_string=='percentile_speed70':
	y_label = '70th percentile of speed \n (Body Length/s)'
if args.a_string=='percentile_speed60':
	y_label = '60th percentile of speed \n (Body Length/s)'
if args.a_string=='percentile_speed100':
	y_label = '100th percentile of speed \n (Body Length/s)'


if args.a_string=='percentile_acc99':
	y_label = '99th percentile of acceleration \n (Body Length/s'+r'$^2$)'

if args.a_string=='percentile_acc90':
	y_label = '90th percentile of acceleration \n (Body Length/s'+r'$^2$)'
if args.a_string=='percentile_acc80':
	y_label = '80th percentile of acceleration \n (Body Length/s'+r'$^2$)'
if args.a_string=='percentile_acc70':
	y_label = '70th percentile of acceleration \n (Body Length/s'+r'$^2$)'
if args.a_string=='percentile_acc60':
	y_label = '60th percentile of acceleration \n (Body Length/s'+r'$^2$)'
if args.a_string=='percentile_acc100':
	y_label = '100th percentile of acceleration \n (Body Length/s'+r'$^2$)'
if args.a_string=='unmasked_startles':
	y_label = 'No. of startles per unit unmasked time (hr' +r'$^{-1}$)'

if args.a_string=='max_loom_speed':
	y_label = 'Maximum speed during loom \n (Body Length/s)'

if args.a_string=='loom_speed99':
	y_label = '99th percentile of speed during loom \n (Body Length/s)'

if args.a_string=='loom_speed90':
	y_label = '90th percentile of speed during loom \n (Body Length/s)'

if args.a_string=='max_loom_acc':
	y_label = 'Maximum acceleration during loom \n (Body Length/s)'+r'$^2$)'

if args.a_string=='loom_acc99':
	y_label = '99th percentile of acceleration during loom \n (Body Length/s)'+r'$^2$)'

if args.a_string=='loom_acc90':
	y_label = '90th percentile of acceleration during loom \n (Body Length/s)'+r'$^2$)'

in_dir1 = '../../output/temp_collective/roi/' + args.a_string + '.p'

annd_values = pickle.load(open(in_dir1, 'rb')) # 'rb is for read binary

in_dir2 = '../../output/temp_collective/roi/' + args.a_string + '_std.p'

out_dir = '../../output/temp_collective/roi_figures/separate_masked_' + args.a_string

std_annd_values = pickle.load(open(in_dir2, 'rb')) # 'rb is for read binary

temperature = [9,13,17,21,25,29]
group = [1,2,4,8,16]
#group = [args.integer_b]
x = 5 # 5 for gs upto 16
#Plotting
lw=1.25
fs=12
colors = plt.cm.viridis(np.linspace(0,1,6))
plt.close('all') # always start by cleaning up
fig = plt.figure(figsize=(8,11))

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=h)
for i in range(6):
    ax = fig.add_subplot(611 + i)
    ax.plot(group[0:x], annd_values[i,0:x], label = str(temperature[i])+ r'$^{\circ}$C', linewidth = lw, color = colors[i])
    ax.fill_between(group[0:x], annd_values[i,0:x] - std_annd_values[i,0:x],  annd_values[i,0:x] + std_annd_values[i,0:x], alpha = 0.3, color = colors[i])
    #plt.xlabel('Group Size', size = fs)
    #plt.ylabel(y_label, size = fs)
    plt.xscale('log',basex=2) 
    if xx == 0:
	    plt.xticks(ticks = [1,2,4,8,16], labels = [1,2,4,8,16])
    else:
	    plt.xticks(ticks = [2,4,8,16], labels = [2,4,8,16])
    #plt.xlim(right = 30)
    ax.tick_params(labelsize=.9*fs)
    #ax.set_title('a)', loc='left', fontsize = fs)
    plt.legend(fontsize=fs, loc='upper right', title = 'Water Temperature', framealpha = 0.5)
fig.add_subplot(111, frame_on=False)

plt.tick_params(labelcolor="none", bottom=False, left=False)


plt.xlabel("Group Size",size = fs*1.5)

plt.ylabel(y_label, size = fs*1.5)

x=6
colors = plt.cm.viridis(np.linspace(0,1,5)) # 5 for gs upto 16
fig1 = plt.figure(figsize=(8,11))

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=h)
for i in range(xx,5):
    if xx == 0:
        ax = fig1.add_subplot(511 + i)
    else:
        ax = fig1.add_subplot(410 + i)
    ax.plot(temperature[0:x], annd_values[0:x,i], label = str(group[i]), linewidth = lw, color = colors[i])
    ax.fill_between(temperature[0:x], annd_values[0:x,i] - std_annd_values[0:x,i],  annd_values[0:x,i] + std_annd_values[0:x,i], alpha = 0.3, color = colors[i])

    #plt.xlabel('Temperature '+r'($^{\circ}$C)', size = fs)
    #plt.locator_params(axis='x', nbins=5)
    #plt.ylabel(y_label, size = fs)
    plt.xticks(ticks = [9,13,17,21,25,29], labels = [9,13,17,21,25,29])
    #plt.xlim(right = 30)
    ax.tick_params(labelsize=.9*fs)
    #ax.set_title('b)', loc='left', fontsize = fs)
    plt.legend(fontsize=fs, loc='upper right', title = 'Group Size', framealpha = 0.5)

fig1.add_subplot(111, frame_on=False)

plt.tick_params(labelcolor="none", bottom=False, left=False)


plt.xlabel('Temperature'+r'($^{\circ}$C)',size = fs*1.5)

plt.ylabel(y_label, size = fs*1.5)

fig.savefig(out_dir+'1.png', dpi = 300)

fig1.savefig(out_dir+'2.png', dpi = 300)


plt.show()
