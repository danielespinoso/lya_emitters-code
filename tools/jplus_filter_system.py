import numpy as np



# map between mock filter-ids and jplus filter-ids
def jpflt(x):
    filt = {
        0 : 'uJAVA',
        1 : 'J0378',
        2 : 'J0395',
        3 : 'J0410',
        4 : 'J0430',
        5 : 'gJAVA',
        6 : 'J0515',
        7 : 'rJAVA',
        8 : 'J0660',
        9 : 'iJAVA',
        10 : 'J0861',
        11 : 'zJAVA',
    }.get(x, ' ')
    if filt == ' ':
        raise ValueError("\n\nfilter band associated to integer: '%s' not found"%x)
        
    return filt




# Dictionary of jplus pivot wavelengths
def jplus_pivot(band='rJAVA'):
    pivot = {
        'uJAVA' : 3522.85761,
        'J0378' : 3786.38436,
        'J0395' : 3950.61134,
        'J0410' : 4100.68634,
        'J0430' : 4300.4497,
        'gJAVA' : 4744.58617,
        'J0515' : 5149.8075,
        'rJAVA' : 6229.77786,
        'J0660' : 6599.77055,
        'iJAVA' : 7676.53339,
        'J0861' : 8602.53784,
        'zJAVA' : 8922.0108022
    }.get(band, ' ')
    if pivot == ' ':
        raise ValueError("\n\nfilter band: '%s' not found"%band)
        
    return pivot
