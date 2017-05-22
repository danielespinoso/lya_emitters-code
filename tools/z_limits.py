
# Dictionary of redshift limits for line emitters
# associated with pass-band waveleghts of jplus
# narrowband filters
def zl(band='J0395', wlength=1215.67):
    passband = {
        # NOTE: at the selected wavelengths include the interval
        # in which the transmission curve is non-zero
        'J0378' : (3649.,3919.),
        'J0395' : (3853.,4028.),
        'J0410' : (3940.,4288.),
        'J0430' : (4135.,4490.),
        'J0515' : (4493.,5298.),
        'J0660' : (6488.,6739.),
        'J0861' : (8289.,8954.),
    }.get(band, ' ')
    if passband[0] == ' ':
        raise ValueError("\n\nfilter band: '%s' not found"%band)
    
    zmin = (passband[0] - wlength)/wlength
    zmax = (passband[1] - wlength)/wlength
    return zmin, zmax
