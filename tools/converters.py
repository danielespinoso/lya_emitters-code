import numpy as np

from .jplus_filter_system import *




# AB magnitudes to fluxes (in ergs) converter
def magtoflux(mags, band='rJAVA', dmags=[], unit='l'):
    c = 2.99792458e18 #in Angstrom/s
    lp = jplus_pivot(band)

    if unit == 'l':
        flux = (10.**(-0.4*(mags + 48.58)))*c/(lp**2.)
    if unit == 'nu':
        flux = 10.**(-0.4*(mags + 48.58))
        
    mask = (mags == -99.)
    flux[mask] = 0.0

    if len(dmags) != 0:
        dflux = flux*(10.**(0.4*dmags)-1.)
        return flux, dflux
    else:
        return flux


        

# fluxes_lambda (in ergs) to AB mags converter
def fluxtomag(flux, band='rJAVA', dflux=[], unit='l'):
    c = 2.99792458e18 #in Angstrom/s
    lp = jplus_pivot(band)
    mags = np.zeros(len(flux))
    
    mask = (flux <= 0)
    mags[mask] = -99.
    
    if unit == 'l':
        mags[~mask] = -2.5*np.log10(((lp**2.)/c)*flux[~mask]) - 48.58
    if unit == 'nu':
        mags[~mask] = -2.5*np.log10(flux[~mask]) - 48.58
    
    if len(dflux) != 0:
        dmags = 2.5*np.log10( (flux+dflux)/flux )
        return mags, dmags
    else:
        return mags




# fluxes_lambda (in ergs) to AB mags converter
def singleMtoF(mags, band='rJAVA', dmags=0, unit='l'):
    c = 2.99792458e18 #in Angstrom/s
    lp = jplus_pivot(band)
    flux = np.empty(1)

    if mags == -99.:
        flux = 0.0
    else:
        if unit == 'l':
            flux = (10.**(-0.4*(mags + 48.58)))*c/(lp**2.)
        if unit == 'nu':
            flux = 10.**(-0.4*(mags + 48.58))
    
    if dmags != 0:
        if mags == -99.:
            dflux = 0.0
        else:
            dflux = flux*(10.**(0.4*dmags)-1.)
        return flux, dflux
    else:
        return flux




# fluxes_lambda (in ergs) to AB mags converter
def singleFtoM(flux, band='rJAVA', dflux=0, unit='l'):
    c = 2.99792458e18 #in Angstrom/s
    lp = jplus_pivot(band)
    mag = np.empty(1)

    if flux <= 0:
        mag = -99.
    else:
        if unit == 'l':
            mag = -2.5*np.log10(((lp**2.)/c)*flux) - 48.58
        if unit == 'nu':
            mag = -2.5*np.log10(flux) - 48.58
    
    if dflux != 0:
        if flux <= 0:
            dmag = -99.
        else:
            dmag = 2.5*np.log10( (flux+dflux)/flux )
        return mag, dmag
    else:
        return mag
