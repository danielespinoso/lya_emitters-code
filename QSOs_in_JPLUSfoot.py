import sys
import deepdish as dd
import numpy as np
import matplotlib.pylab as plt

import tools
setup = tools.set_up()

sys.path.append(setup['jplus_code'])
import jplus
from jplus.tools import crossmatch_angular as xmatch


#---------------- LOAD DATSETS ----------------#
#qso = jplus.datasets.fetch_sdss_qso(mag_type="aperMags", overwrite=False, mag_limit=[14,23])
qso = dd.io.load(setup['sdss_qso'])     # loads sdss qso catalogue
jpl = dd.io.load(setup['jplus_input'])  # loads the initial jplus catalogue
gaia = dd.io.load(setup['gaia'])        # loads gaia catalogue


#------------------- Z-MASK -------------------#
zmin, zmax = tools.zl(setup['filters'][0])
zmask_qso = ( (qso['zspec'] > zmin) & (qso['zspec'] < zmax))
print '\nNumber of all right-z QSOs:  ', len(qso['coords'][zmask_qso])


#------- PLOTS JPLUS FOOTPRINT and COUNT QSOs in IT --------#
qsos_in_jpfoot = [],[]
fig = plt.figure(figsize=(12,10))
for i in set(jpl['tile_id']):
    
    tilemask = (jpl['tile_id'] == i)
    minra = min(jpl['coords'][tilemask,0]); maxra = max(jpl['coords'][tilemask,0])
    mindec = min(jpl['coords'][tilemask,1]); maxdec = max(jpl['coords'][tilemask,1])

    coordmask = ((qso['coords'][:,0] > minra) & (qso['coords'][:,0] < maxra) & \
                 (qso['coords'][:,1] > mindec) & (qso['coords'][:,1] < maxdec))
    totmask = ((coordmask) & (zmask_qso))

    #select QSOs in JPLUS footprint excluding doubles in overlap regions between tiles
    for x, y in zip(qso['coords'][totmask,0], qso['coords'][totmask,1]):
        if x not in qsos_in_jpfoot[0] and y not in qsos_in_jpfoot[1]:
            qsos_in_jpfoot[0].append(x)
            qsos_in_jpfoot[1].append(y)
    #plot jplus footprint (one square per tile)
    plt.plot( (minra, maxra), (mindec, mindec), 'k-', alpha=0.6 )
    plt.plot( (maxra, maxra), (mindec, maxdec), 'k-', alpha=0.6 )
    plt.plot( (maxra, minra), (maxdec, maxdec), 'k-', alpha=0.6 )
    plt.plot( (minra, minra), (maxdec, mindec), 'k-', alpha=0.6 )
#plt.plot(qsos_in_jpfoot[0], qsos_in_jpfoot[1], 'or', markersize=4, alpha=0.3)

nQz_JPfoot = len(qsos_in_jpfoot[0])  #number of right-z SDSS QSOs in JPLUS footprint
print 'Number of right-z QSOs in JPLUS footprint:  ', nQz_JPfoot


#-------   JPLUS - SDSS_QSOs(right_z) CROSSMATCH   --------#
dist, ind = xmatch(jpl['coords'], qso['coords'][zmask_qso], max_distance=3./3600.)
qso_mask = (dist != np.inf)


retrived = len(jpl['coords'][qso_mask,0])
print 'Number of RETRIVED right-z QSOs in JPLUS catalogue:  ', retrived
print 'Fraction of RETRIVED right-z QSOs:  ', int(10000*float(retrived)/nQz_JPfoot)/100., '%'

if setup['plot_gaia'] == True:
    plt.plot(gaia['coords'][:,0], gaia['coords'][:,1], 'oc', markersize=1, alpha=0.3)
#plt.plot(jpl['coords'][qso_mask,0],jpl['coords'][qso_mask,1], 'ob', markersize=2, alpha=0.3)
plt.show()
plt.close()
sys.exit()
