import numpy as np
import sys
import os
import deepdish as dd

import tools
setup = tools.set_up()
import tools.screen_info as tsi

sys.path.append(setup['jplus_code'])
import jplus
from jplus.datasets.sdss_objects import *

def get_sdss_galaxies(num_of_chuncks=20):
    nch = num_of_chuncks
    if not os.path.exists(setup['sdss_gal']):
        tsi.update_print('Downloading SDSS galaxies... ', appendix=' ')
        sd_gal = fetch_sdss_objects(mag_type="modelMags", overwrite=False, object_name="galaxies", casjobs=True, mag_limit=[0,21], spectroscopic=True, nchunks=nch)
        tsi.update_print('Downloading SDSS galaxies... ', appendix='done')
    else:
        tsi.update_print('loading SDSS galaxies... ', appendix=' ')
        sd_gal = dd.io.load(setup['sdss_gal'])
        tsi.update_print('loading SDSS galaxies... ', appendix='done')
    return sd_gal



def get_sdss_qso(num_of_chuncks=20):
    nch = num_of_chuncks
    if not os.path.exists(setup['sdss_qso']):
        tsi.update_print('Downloading SDSS quasars... ', appendix=' ')
        sd_qso = fetch_sdss_objects(mag_type="modelMags", overwrite=False, object_name="qso", casjobs=True, mag_limit=[14,23], spectroscopic=True, nchunks=nch)
        tsi.update_print('Downloading SDSS quasars... ', appendix='done')
    else:
        tsi.update_print('loading SDSS quasars... ', appendix=' ')
        sd_qso = dd.io.load(setup['sdss_qso'])
        tsi.update_print('loading SDSS quasars... ', appendix='done')
    return sd_qso



def get_sdss_stars(num_of_chuncks=20):
    nch = num_of_chuncks
    if not os.path.exists(setup['sdss_str']):
        tsi.update_print('Downloading SDSS stars... ', appendix=' ')
        sd_stars = fetch_sdss_objects(mag_type="modelMags", overwrite=False, object_name="stars", casjobs=True, mag_limit=[10,17], spectroscopic=True, nchunks=nch)
        tsi.update_print('Downloading SDSS stars... ', appendix='done')
    else:
        tsi.update_print('loading SDSS stars... ', appendix=' ')
        sd_stars = dd.io.load(setup['sdss_str'])
        tsi.update_print('loading SDSS stars... ', appendix='done')
    return sd_stars
    
