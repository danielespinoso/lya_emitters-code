import sys
import time
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt


#####################  ATTENTION  ####################
# This code reads a general CSV input and formats
# it into a python dictionary (each column becomes a
# key in the dictionary). Obviously the exact structure
# of the input file can differ each time, thus some
# parameters/variables must be enabled/disabled
# manually depending on the input.
######################################################


#----------------------#
#       SETTINGS       #
#----------------------#
inp = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/gaia_backups/gaia_alljplfootage_(patched)_new.csv'
out = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/gaia_jpFoot_new.h5'


#------------------------#
#    LOAD AND CONVERT    #
#------------------------#
print 'reading output...'
data = np.genfromtxt(inp, delimiter=',', skip_header=1, unpack=False) # reading actual data
print 'output read'

print 'saving output...'
table = {
    'coords' : np.array([data[:,1],data[:,2]]).T
}

dd.io.save(out, table)
