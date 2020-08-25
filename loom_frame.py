"""
Goal - to figure out which videos do not have the loom frame number associated with the group size and replicate
"""
import csv
import numpy as np

temperature = range(9,30,4)

group = [1,2,4,8,16]

replication = range(10) # number of replicates per treatment

ii = 0 # to keep count of temperature



rows = []
with open('../../data/temp_collective/looms_roi.csv', 'r') as csvfile:
    looms = csv.reader(csvfile)
    for row in looms:
        rows.append(row)
        

for temp in temperature:
    jj = 0 # to keep count of groups
    for j in group:
    	for k in replication:
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
    		g = str(j)
    		r = str(k+1)
    		loom = np.zeros([5,1])        
    		for i in range(len(rows)):
    		    if rows[i][1]==cam and rows[i][3]==g and rows[i][4]==r:
    		        for m in range(5):
    		            loom[m] = int(rows[i][2]) + m*11403 
    		if loom[0] == 0:
    			print(temp,j,k+1)
