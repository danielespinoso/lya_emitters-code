import numpy as np

def linfit(x, y, sigma):
    w = 1./(sigma**2.)
    Sw = sum(w)
    
    wx = x/(sigma**2.)
    Swx = sum(wx)
    
    wy = y/(sigma**2.)
    Swy = sum(wy)
    
    wxy = (x*y)/(sigma**2.)
    Swxy = sum(wxy)
    
    x2 = x**2.
    wx2 = x2/(sigma**2.)
    Swx2 = sum(wx2)

    D = ((Sw) * (Swx2)) - ((Swx)**2.)

    nA = ((Swx2)*(Swy)) - ( (Swx)*(Swxy) )
    nB = (Sw)*(Swxy) - ( (Swx)*(Swy) )
    A = nA/D
    B = nB/D
    
    return A,B
