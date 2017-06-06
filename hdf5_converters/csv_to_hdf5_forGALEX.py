import sys
import time
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt


#################  ATTENTION  #################
# This code reads and formats galex csv catalogue. It also
# cuts it by keeping only sources in jplus footprint
# Be careful with the numbers of columns
# of the input data (they depend on the
# exact query submitted to the jplus
# archive). The 'table - data' fields
# (below) must be adapted consequently!!
###############################################

jplus_input = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/jplus_allTOr24_autoMags_upad_dual_jplusPhot.h5'
alljpl = dd.io.load(jplus_input)

# input galex filename
infn = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/galex_photoObjs_autoMags.csv'

# output galex filename
archive_file = '/home/CEFCA/dspinoso/works/lya_emitters/datasets/galex_photoObjs_autoMags.h5'


# query text (the exact text is optional, at least put a non-empty string)
query_text = 'SELECT\nra, dec, NUV_MAG_AUTO, NUV_MAGERR_AUTO, FUV_MAG_AUTO, FUV_MAGERR_AUTO\ninto mydb.MyTable_1\nfrom\nGALEX_GR6plus7.photoobjall\nWHERE\nnuv_s2n > 3 AND fuv_s2n > 3'

print 'reading data...'
data = np.genfromtxt(infn, delimiter=',', skip_header=1, unpack=False ) 

table = {'coords':data[:,[0,1]],
         'nuv':data[:,[2,3]],
         'fuv':data[:,[4,5]],
         'SQL_query':query_text,
         'date':time.strftime("%c"),
         'filename':archive_file
}


# SELECT ONLY GALEX SOURCES IN JPLUS FOOTPRINT
dummy = [],[]
nummy = [],[]
fummy = [],[]
for i in set(alljpl['tile_id']):
    tilemask = (alljpl['tile_id'] == i)
    minra = min(alljpl['coords'][tilemask,0]); maxra = max(alljpl['coords'][tilemask,0])
    mindec = min(alljpl['coords'][tilemask,1]); maxdec = max(alljpl['coords'][tilemask,1])
    coordmask = ((table['coords'][:,0] > minra) & (table['coords'][:,0] < maxra) & \
                 (table['coords'][:,1] > mindec) & (table['coords'][:,1] < maxdec))
    dummy[0].extend(table['coords'][coordmask,0])
    dummy[1].extend(table['coords'][coordmask,1])
    nummy[0].extend(table['nuv'][coordmask,0])
    nummy[1].extend(table['nuv'][coordmask,1])
    fummy[0].extend(table['fuv'][coordmask,0])
    fummy[1].extend(table['fuv'][coordmask,1])
    #plt.plot(alljpl['coords'][tilemask,0],alljpl['coords'][tilemask,1], 'ob', markersize=2, alpha=0.3)

table['coords'] = np.array(dummy).T
table['nuv'] = np.array(nummy).T
table['fuv'] = np.array(fummy).T

#print  table['coords'].shape
#plt.plot(table['coords'][:,0],table['coords'][:,1], 'or', markersize=4, alpha=0.3)
#plt.show()
#plt.close()
#sys.exit()

print '\nsaving to hdf5 file...'

dd.io.save(archive_file, table)

print '\ndone'
