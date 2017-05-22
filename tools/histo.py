
import numpy as np


####  HISTOGRAM (by David Izquierdo) ####
# makes an histogram with x-values in the center of the x-bins
# (improvement of numpy.histogram)
def histogram(x,xmin,xmax,sizebins):
    bins=np.arange(xmin,xmax,sizebins)
    hist, b_edges = np.histogram(x,bins)
    centers = (b_edges[:-1] + b_edges[1:]) / 2
    hist = np.array(hist, dtype=np.float32)
    return centers, hist
