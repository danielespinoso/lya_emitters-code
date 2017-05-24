import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import tools.settings
setup = tools.settings.set_up()

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gauss_double(x, *p):
    A, mu, sigma, B, mv, sogma = p
    return (A*np.exp(-(x-mu)**2/(2.*sigma**2))) + (B*np.exp(-(x-mv)**2/(2.*sogma**2)))

def magerr_fit(x, *p):
    a, b, c = p
    return np.sqrt( a*(10.**(b*x+c)) )



# THE FOLLOWING FUNCTION TRIES TO FIT A SINGLE OR DOUBLE GAUSSIAN TO THE
# "MORPHOLOGY CLOUD" in the [(mu_max - rJAVA)  vs.  (rJAVA)] plane.
# It proceeds on magnitude bins ("slices"): the first attempt is always
# fitting a single gaussian to the (mu_max - rJAVA) histogram in that bin.
# When the fitted gaussian is too broad, the histogram get splitted and
# the code tries to fit a double gaussian. As a last sep, if the two fitted
# gaussians overlap, then the first "broad" gaussian is retrived as "best fit"

def evaluate_morphology(data, extent, mask_tile, tilenum):
    lims = [15., 21.5]
    binlen = 0.2
    bb = 80
    fact = 2.
    plot_by_slice = False

    magbins = []      #mag vector to plot first gaussian fits
    magbins_2 = []    #mag vector to plot second gaussian fits (where defined)
    mean_magerr = []  #mean error on magnitudes as a function of magnitudes
    amp_pb_1 = []     #amplitude of the first gaussian fit
    mean_pb_1 = []    #first-gaussian central-value at each mag_bin
    stdv_pb_1 = []    #standard deviation of the first gaussian distribution
    mean_pb_2 = []    #second-gaussian central-value at each mag_bin
    stdv_pb_2 = []    #standard deviation of the second gaussian distribution

    mm = ((data['rJAVA'][:,0] > 15.) & (data['rJAVA'][:,0] < 16.) & (mask_tile))
    fr = np.mean(extent[mm])  #first reference (tile dependent)
    ylims = (fr-1.5, fr+3.5)
    
    for i in np.arange(lims[0], lims[1], binlen):
        sliced = ((data['rJAVA'][:,0] > i) & (data['rJAVA'][:,0] < i+binlen) & (mask_tile))
        if len(extent[sliced]) > 30:
            
            hist, bin_edges = np.histogram(extent[sliced], range=ylims, bins=bb)
            bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
            avg = np.mean(extent[sliced])
            
            p0 = [5., avg, 0.2]   # initial guess on the gaussian-fit parameters
            coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0) # gaussian fit
            amp = coeff[0]        # get the amplitude
            mean = coeff[1]       # get the mean
            stdv = abs(coeff[2])  # get the st.dev
            
            if stdv/amp < 0.09: #check on the width of the first gaussian-fit
                mean_pb_1.append(mean)
                stdv_pb_1.append(stdv)
                
                if plot_by_slice == True:
                    hist_fit = gauss(bin_centres, *coeff)   #only to plot the gaussian fit
                    plt.plot(bin_centres, hist_fit, 'r--', label='Gauss fit')
                    plt.hist(extent[sliced], range=ylims, bins=bb, color='b', alpha=0.5, label='Data')
                    plt.show()
                    plt.close()
                
            else:  #if the first gaussian-fit is too wide then try with two gaussians
                upto = mean_pb_1[0]+0.5
                hist_1, bin_edges_1 =np.histogram(extent[sliced],range=(ylims[0], upto),bins=bb/2)
                bin_centres_1 = (bin_edges_1[:-1] + bin_edges_1[1:])/2
                
                p0_1 = [10., 1.8, 0.2]
                coeff_1, var_matrix = curve_fit(gauss, bin_centres_1, hist_1, p0=p0_1) #1 gauss fit
                mean_1 = coeff_1[1]
                stdv_1 = abs(coeff_1[2])

                #-- SECOND GAUSSIAN --#
                hist_2, bin_edges_2 =np.histogram(extent[sliced],range=(upto, ylims[1]),bins=bb/2)
                bin_centres_2 = (bin_edges_2[:-1] + bin_edges_2[1:])/2
                
                p0_2 = [5., 3., 0.2]
                coeff_2, var_matrix = curve_fit(gauss, bin_centres_2, hist_2, p0=p0_2) #fit
                mean_2 = coeff_2[1]
                stdv_2 = abs(coeff_2[2])

                #check the separation of the two gaussians
                separ = abs(mean_2 - mean_1)
                check = stdv_1# + stdv_2
                if separ < check:
                    mean_pb_1.append(mean) #retrives the original mean (first single gaussian fit)
                    stdv_pb_1.append(stdv) #retrives the original stdv (first single gaussian fit)
                    separ_limit = i
                    
                    if plot_by_slice == True:
                        # #---- PLOT BY SLICE ----#
                        coeff = [amp, mean, stdv]
                        hist_fit = gauss(bin_centres, *coeff)   #only to plot the gaussian fit
                        plt.plot(bin_centres, hist_fit, 'r--', label='Gauss fit')

                else:
                    mean_pb_1.append(mean_1)
                    stdv_pb_1.append(stdv_1)
                    mean_pb_2.append(mean_2)
                    stdv_pb_2.append(stdv_2)
                    magbins_2.append(i)
                    
                    if plot_by_slice == True:
                        # #---- PLOT BY SLICE ----#
                        hist_fit_1 = gauss(bin_centres_1, *coeff_1)   #only to plot the gaussian fit
                        plt.plot(bin_centres_1, hist_fit_1, 'r--', label='Gauss fit')
                        hist_fit_2 = gauss(bin_centres_2, *coeff_2)   #only to plot the gaussian fit
                        plt.plot(bin_centres_2, hist_fit_2, 'c--', label='Gauss fit')
                
                        plt.hist(extent[sliced], range=ylims, bins=bb, color='b', alpha=0.5)
                        plt.show()
                        plt.close()
            
            magbins.append(i)
            amp_pb_1.append(amp)
            mean_magerr.append(np.mean(data['rJAVA'][sliced,1]))
            
    magbins = np.array(magbins)
    magbins_2 = np.array(magbins_2)
    mean_pb_1 = np.array(mean_pb_1)
    stdv_pb_1 = np.array(stdv_pb_1)
    mean_pb_2 = np.array(mean_pb_2)
    stdv_pb_2 = np.array(stdv_pb_2)
    mean_magerr = np.array(mean_magerr)

    amp_pb_1 = np.array(amp_pb_1)

    # p1 = [1., 0.4, 40.]
    # pars, varmatr = curve_fit(magerr_fit, magbins, mean_magerr, p0=p1)
    # fitted_magerr = magerr_fit(mean_magerr, *pars)
    
    tot_err = np.empty(len(mean_magerr))
    for i, j in zip(mean_magerr, range(len(mean_magerr))):
        if stdv_pb_1[j] < 0.3:
            tot_err[j] = np.sqrt(mean_magerr[j]**2. + stdv_pb_1[j]**2.)
        else:
            tot_err[j] = mean_magerr[j]
    
    magbin_mask = (magbins <= 18.0)
    a, b = tools.linfit(magbins[magbin_mask], mean_pb_1[magbin_mask], stdv_pb_1[magbin_mask])
    #ix = np.arange(12., 22., 0.01)
    #iy = a+b*ix
    iy = a+b*magbins

    plt.plot(magbins, stdv_pb_1/amp_pb_1)
    plt.title('Tile: ' + str(int(tilenum)))
    plt.savefig(setup['plots']+'stdev_over_amplitude_tile'+str(tilenum)+'.png')
    plt.close()

    f = plt.figure(figsize=(12,10))
    plt.plot(data['rJAVA'][mask_tile,0], extent[mask_tile], 'og', markersize=5, alpha=0.1)
    plt.errorbar(magbins, mean_pb_1, yerr=stdv_pb_1, fmt='o', ecolor='r', color='r', markersize=2, linewidth=2)
    plt.errorbar(magbins_2, mean_pb_2, yerr=stdv_pb_2, fmt='o', ecolor='m', color='m', markersize=2, linewidth=2)
    plt.plot(magbins, iy, 'k-')
    plt.plot(magbins, iy+mean_magerr, 'b-')
    plt.plot(magbins, iy-mean_magerr, 'b-')
    plt.plot(magbins, iy+tot_err, 'c-')
    plt.plot(magbins, iy-tot_err, 'c-')
    plt.title('Tile: ' + str(int(tilenum)))
    plt.xlabel('rJAVA  [mags]')
    plt.ylabel('mu_max - rJAAVA  [auto-mags]')
    plt.savefig(setup['plots']+'interesting_morpho-plot_tile'+str(tilenum)+'.png')
    #plt.show()
    plt.close()
    









#------- DOWN HERE: AN ALMOST-WORKING ATTEMPT OF DOUBLE-GAUSSIAN FIT -------#

# def evaluate_morphology(data, extent, mask_tile):
#     lims = [15., 22.]
#     binlen = 0.2
#     fact = 2.

#     magbins_1 = []
#     magbins_2 = []
#     mean_pb_1 = []  #mean per mag_bin of the first gaussian (extendedness cloud gaussian fit)
#     mean_pb_2 = []  #mean per mag_bin of the second gaussian (extendedness cloud gaussian fit)
#     stdv_pb_1 = []  #standard deviation of the first gaussian
#     stdv_pb_2 = []  #standard deviation of the second gaussian

#     mm = ((data['rJAVA'][:,0] > 15.) & (data['rJAVA'][:,0] < 16.) & (mask_tile))
#     fr = np.mean(extent[mm])  #first reference (tile dependent)
#     ylims = (fr-1., fr+2)
    
#     for i in np.arange(lims[0], lims[1], binlen):
#         sliced = ((data['rJAVA'][:,0] > i) & (data['rJAVA'][:,0] < i+binlen) & (mask_tile))
#         if len(extent[sliced]) > 50:
        
#             hist, bin_edges = np.histogram(extent[sliced], range=ylims, bins=100)
#             bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
#             avg = np.mean(extent[sliced])

#             if i < 19.0 or i >= 20.5:
#                 p0 = [5., avg, 0.2]   # initial guess on the gaussian-fit parameters
#                 coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0) # gaussian fit
#                 mean = coeff[1]       # get the mean
#                 stdv = abs(coeff[2])  # get the st.dev

#                 mean_pb_1.append(coeff[1])
#                 stdv_pb_1.append(coeff[2])
                            
#                 # # #---- PLOT BY SLICE ----#
#                 # hist_fit = gauss(bin_centres, *coeff)   #only to plot the gaussian fit
#                 # plt.hist(extent[sliced], range=ylims, bins=100, color='b', alpha=0.5, label='Data')
#                 # plt.plot(bin_centres, hist_fit, 'r--', label='Gauss fit')
#                 # plt.show()
#                 # plt.close()

#             else:
#                 p0 = [5., avg, 0.2, 1., (avg+2), 0.5]   # initial guess on the gaussian-fit parameters
#                 coeff, var_matrix = curve_fit(gauss_double, bin_centres, hist, p0=p0) # gaussian fit
#                 mean = coeff[1]       # get the mean
#                 stdv = abs(coeff[2])  # get the st.dev
                
#                 mean_pb_1.append(coeff[1])
#                 stdv_pb_1.append(coeff[2])
#                 if coeff[4] > 1. and coeff[4] < 4.:
#                     magbins_2.append(i)
#                     mean_pb_2.append(coeff[4])
#                     stdv_pb_2.append(coeff[5])
            
#                 # # #---- PLOT BY SLICE ----#
#                 # hist_fit = gauss_double(bin_centres, *coeff)   #only to plot the gaussian fit
#                 # plt.hist(extent[sliced], range=ylims, bins=100, color='b', alpha=0.5, label='Data')
#                 # plt.plot(bin_centres, hist_fit, 'r--', label='Gauss fit')
#                 # plt.show()
#                 # plt.close()

#             magbins_1.append(i)

#     magbins_1 = np.array(magbins_1)
#     magbins_2 = np.array(magbins_2)
#     mean_pb_1 = np.array(mean_pb_1)
#     stdv_pb_1 = np.array(stdv_pb_1)
#     mean_pb_2 = np.array(mean_pb_2)
#     stdv_pb_2 = np.array(stdv_pb_2)

#     a, b = tools.linfit(magbins_1, mean_pb_1, stdv_pb_1)
#     c, d = tools.linfit(magbins_2, mean_pb_2, stdv_pb_2)
#     ix = np.arange(12., 24., 0.01)
#     iy = a+b*ix
#     iz = c+d*ix
    
#     f = plt.figure(figsize=(12,10))
#     plt.plot(data['rJAVA'][mask_tile,0], extent[mask_tile], 'og', markersize=5, alpha=0.1)
#     plt.errorbar(magbins_1, mean_pb_1, yerr=stdv_pb_1, fmt='o', ecolor='r', color='r', markersize=2, linewidth=2)
#     plt.errorbar(magbins_2, mean_pb_2, yerr=stdv_pb_2, fmt='o', ecolor='b', color='b', markersize=2, linewidth=2)
#     plt.plot(ix, iy, 'k-')
#     plt.plot(ix, iz, 'b-')
#     plt.show()
#     plt.close()
