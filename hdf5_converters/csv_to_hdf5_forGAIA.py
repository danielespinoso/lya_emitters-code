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
inp = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/gaia_backups/gaia_alljplfootage_(patched).csv'
out = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/gaia_jpFoot_new.h5'


#------------------------#
#    LOAD AND CONVERT    #
#------------------------#
print 'reading input...'
data = np.genfromtxt(inp, delimiter=',', skip_header=1, unpack=False) # reading actual data
print 'input read'

query_text = 'SELECT\n`ra`, `dec`\nFROM\nGDR1.gaia_source\nWHERE\n`ra` BETWEEN 126.25 AND 142.42 AND\n`dec` BETWEEN 51 AND 60 OR\n`ra` BETWEEN 179.70 AND 245.30 AND\n`dec` BETWEEN 53.35 AND 57.6 OR\n`ra` BETWEEN 257.3 AND 267.6 AND\n`dec` BETWEEN 56.1 AND 57.6 OR\n`ra` BETWEEN 12.25 AND 16.5 AND\n`dec` BETWEEN 38.6 AND 41.6 OR\n`ra` BETWEEN 21.3 AND 23.7 AND\n`dec` BETWEEN 38.7 AND 41.6 OR\n`ra` BETWEEN 0.0 AND 3.58 AND\n`dec` BETWEEN 35.9 AND 38.8 OR\n`ra` BETWEEN 21.4 AND 28.2 AND\n`dec` BETWEEN 33.16 AND 34.6 OR\n`ra` BETWEEN 33.7 AND 40.0 AND\n`dec` BETWEEN 29.00 AND 33.24 OR\n`ra` BETWEEN 104.6 AND 171.9 AND\n`dec` BETWEEN 36.62 AND 42.30 OR\n`ra` BETWEEN 114.8 AND 116.9 AND\n`dec` BETWEEN 45.0 AND 46.5 OR\n`ra` BETWEEN 107.90 AND 153.90 AND\n`dec` BETWEEN 29.7 AND 32.55 OR\n`ra` BETWEEN 247.25 AND 249.17 AND\n`dec` BETWEEN 39.45 AND 40.91 OR\n`ra` BETWEEN 272.22 AND 285.29 AND\n`dec` BETWEEN 38.04 AND 42.33 OR\n`ra` BETWEEN 257.60 AND 282.85 AND\n`dec` BETWEEN 29.69 AND 34.00 OR\n`ra` BETWEEN 254.2 AND 259.15 AND\n`dec` BETWEEN 22.72 AND 27.02 OR\n`ra` BETWEEN 329.95 AND 336.51 AND\n`dec` BETWEEN 26.15 AND 31.85 OR\n`ra` BETWEEN 334.85 AND 348.38 AND\n`dec` BETWEEN 31.75 AND 34.63 OR\n`ra` BETWEEN 344.96 AND 349.57 AND\n`dec` BETWEEN 22.03 AND 23.50 OR\n`ra` BETWEEN 18.80 AND 39.00 AND\n`dec` BETWEEN 42.8 AND 45.76 OR\n`ra` BETWEEN 104.40 AND 119.67 AND\n`dec` BETWEEN 72.85 AND 79.93\n'

print 'saving output...'
table = {
    'query' : query_text,
    'coords' : np.array([data[:,1],data[:,2]]).T
}

dd.io.save(out, table)
