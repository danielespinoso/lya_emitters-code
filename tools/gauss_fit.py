
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gauss_double(x, *p):
    A, mu, sigma, B, mv, sogma = p
    return (A*np.exp(-(x-mu)**2/(2.*sigma**2))) + (B*np.exp(-(x-mv)**2/(2.*sogma**2)))


# This function fits a gaussian to a histrogram.
# Once the data are given it produces the histogram
# according to the range and number of bins then fits it
# with a gaussian starting from the initial guess for
# amplitude, central value and sigma (the 'init_guess' input
# must be ordered this way: ampl., mean, sigma)
# RETURNS AMPLITUDE, MEAN AND SIGMA OF THE FIT (in this order)
def fit_gauss_to_histogram(data, hist_range, bin_num, init_guess, guess_from_avg=True, show_plot=True):

    hist, bin_edges = np.histogram(data, range=hist_range, bins=bin_num)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    if guess_from_avg == True:
        avg = np.mean(data)
        init_guess[1] = avg
    
    p0 = init_guess       # initial guess on the gaussian-fit parameters
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0) # gaussian fit
    
    A = coeff[0]       # get the amplitude
    B = coeff[1]       # get the mean
    C = abs(coeff[2])  # get the st.dev
            
    if show_plot == True:
        hist_fit = gauss(bin_centres, *coeff)   #only to plot the gaussian fit
        plt.plot(bin_centres, hist_fit, 'r--', label='Gauss fit')
        plt.hist(data, range=hist_range, bins=bin_num, color='b', alpha=0.5)
        plt.legend()
        plt.show()
        plt.close()
    
    return A, B, C
