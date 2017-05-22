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
inp = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/jplus_tile_image_(03-05-2017).csv'
out = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/jplus_tile_image_(03-05-2017).h5'


#------------------------#      Only for .csv data retrived from JPLUS catalogue:
#    LOAD AND CONVERT    #         - the first row contains the query text used to retrive the data
#------------------------#         - the second row contains the names of the columns

with open(inp, 'r') as f:
    query_text = f.readline()  # reading input first line
    names = f.readline()       # reading input second line
f.close()

names = names.lower()
kys = names[1:-2].split('","') # formatting dictionary keys

data = np.genfromtxt(inp, delimiter=',', skip_header=2, unpack=False) # reading actual data

table = {}
for i in range(len(kys)):
    table[kys[i]] = data[:,i]  # writing the dictionary
table['tile_idd'] = table.pop('tile_id')

from tools.jplus_filter_system import jpflt
table['filt_id'] = []
for i in table['filter_id']:   # converting filter indexes to filter names
    table['filt_id'].append( jpflt( int(i - 1)) )
table['filt_id'] = np.array(table['filt_id']).T

# changing keys names to more intuitive ones
table['psf'] = table.pop('fwhm_mean')
table['fwhm_gauss'] = table.pop('fwhm')
table['fwhm_gauss_rms'] = table.pop('fwhmg_rms')
table['tile_id'] = table.pop('ref_tile_id')
table['image_id'] = table.pop('tile_idd')
table['filter_index'] = table.pop('filter_id')
table['filter_name'] = table.pop('filt_id')

dd.io.save(out, table)
