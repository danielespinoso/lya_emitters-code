import os
import sys
import numpy as np
import pyfits as pf

import tools.settings
setup = tools.settings.set_up()

def read_sdss_spectrum(speck):
    lines = pf.open(speck)
    lines.close()

    print lines[0].header
    sys.exit()
    c0 = lines[0].header['coeff0']
    c1 = lines[0].header['coeff1']
    npix = lines[1].header['naxis2']
    units = lines[0].header['bunit']
    
    lamb = 10.**(c0 + c1 * np.arange(npix))
    flux = lines[1].data

    return lamb, flux, units
