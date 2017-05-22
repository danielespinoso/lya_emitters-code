import sys
import time
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt


#####################  ATTENTION  ####################
# This code reads a general ascii input (a text table)
# and formats it into a python dictionary (column are
# rearranged into dictionary keys. Obviously the exact
# Structure of the input file can differ each time,
# thus some parameters/variables must be enabled/disabled
# manually depending on the input.
######################################################


#----------------------#
#       SETTINGS       #
#----------------------#
inp = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/daniele.dat'
out = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/rafa_HII.h5'


#------------------------#
#    LOAD AND CONVERT    #
#------------------------#

data = np.genfromtxt(inp, skip_header=1, usecols=(0,1,2), unpack=False) # reading actual data

for i in range(len(data[:,0])):
    if data[i,0].dtype != np.float64:
        print 'NO'
        
table = {
    'coords' : np.array([data[:,0],data[:,1]]).T,
    'pixels_over_threshold' : np.array(data[:,2]).T
}

dd.io.save(out, table)
