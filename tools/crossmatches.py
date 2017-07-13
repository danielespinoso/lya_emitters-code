import os
import sys
import numpy as np
import deepdish as dd

import tools.settings
setup = tools.settings.set_up()
from tools.screen_info import update_print as uPrint
    
sys.path.append(setup['jplus_code'])
import jplus
from jplus.tools import crossmatch_angular as ang_xmatch

def xmatch_wrapper(jpldata):
    
    if setup['plot_sdssGal'] == True:
        database = tools.get_sdss_galaxies()  #loads and crossmatch sdss galaxies catalogue
        uPrint('crossmatching SDSS galaxies....      ', appendix='')
        dis, ind = ang_xmatch(jpldata['coords'], database['coords'], max_distance=3/3600.)
        galaxies = (dis != np.inf)
        uPrint('crossmatching SDSS galaxies....      ', appendix='done')
    else:
        galaxies = np.ones(len(jpldata['coords'][:,0]), dtype=bool)
    

    if setup['plot_sdssQSO'] == True:
        database = tools.get_sdss_qso()       #loads sdss quasars catalogue
        uPrint('crossmatching SDSS quasars....       ', appendix='')
        dis, ind = ang_xmatch(jpldata['coords'], database['coords'], max_distance=3/3600.)
        quasars = (dis != np.inf)
        
        zetamask = ((database['zspec'] > setup['zmin']) & (database['zspec'] < setup['zmax']))
        database_z = jplus.tools.select_object(database, zetamask)
        dis, ind = ang_xmatch(jpldata['coords'], database_z['coords'], max_distance=3/3600.)
        quasars_z = (dis != np.inf)
        uPrint('crossmatching SDSS quasars....       ', appendix='done')
    else:
        quasars = np.ones(len(jpldata['coords'][:,0]), dtype=bool)
        quasars_z = np.ones(len(jpldata['coords'][:,0]), dtype=bool)

        
    if setup['plot_sdssSTAR'] == True:
        database = tools.get_sdss_stars()     #loads sdss stars catalogue
        uPrint('crossmatching SDSS stars....         ', appendix='')
        dis, ind = ang_xmatch(jpldata['coords'], database['coords'], max_distance=3/3600.)
        stars = (dis != np.inf)
        uPrint('crossmatching SDSS stars....         ', appendix='done')
    else:
        stars = np.ones(len(jpldata['coords'][:,0]), dtype=bool)
        
    
    if setup['plot_gaia'] == True:
        uPrint('loading gaia... ', appendix=' ')
        database = dd.io.load(setup['gaia'])      #loads gaia stars catalogue
        uPrint('loading gaia....      ', appendix='done')
        uPrint('crossmatching gaia....      ', appendix='')
        dis, ind = ang_xmatch(jpldata['coords'], database['coords'], max_distance=3/3600.)
        gaia = (dis != np.inf)
        uPrint('crossmatching gaia....      ', appendix='done')
    else:
        gaia = np.ones(len(jpldata['coords'][:,0]), dtype=bool)

    return galaxies, quasars, quasars_z, stars, gaia
