import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import tools.settings
setup = tools.settings.set_up()
from tools.gauss_fit import fit_gauss_to_histogram as FGtH

def magerr_fit(x, *p):
    a, b, c = p
    return a + (10.**(b*x+c))


# THE FOLLOWING FUNCTION TRIES TO FIT A SINGLE GAUSSIAN TO THE
# "POINTLIKE MORPHOLOGY CLOUD" in the [(mu_max - rJAVA)  vs.  (rJAVA)] plane.
# It proceeds on magnitude bins ("slices"): as a first attempt it always
# tries to fit a single gaussian to the (mu_max - rJAVA) histogram in that bin.
# When the fitted gaussian is too broad, the histogram get splitted and
# the code tries to fit a gaussian only on the "pointlike" cloud. When the two
# clouds are indistinguishable the first fit is usually tight enough to avoid
# the splitting of the histogram so the pointlike and extended clouds get fitted
# as a single distribution

def evaluate_morphology(data, extent, mask_tile, tilenum):
    lims = [16., 20.5]
    binlen = 0.2
    bb = 80
    fact = 2.
    plot_by_slice = False

    magbins = []      #mag vector to plot first gaussian fits
    mean_magerr = []  #mean error on magnitudes as a function of magnitudes
    amp_pb = []       #amplitude of the first gaussian fit
    mean_pb = []      #first-gaussian central-value at each mag_bin
    stdv_pb = []      #standard deviation of the first gaussian distribution

    mm = ((data['mag_auto_r'][:] > 16.) & (data['mag_auto_r'][:] < 16.+binlen) & (mask_tile))
    fr = np.mean(extent[mm])  #first reference (tile dependent)
    ylims = (fr-2.5, fr+3.5)
    
    for i in np.arange(lims[0], lims[1], binlen):
        sliced = ((data['mag_auto_r'][:] > i) & (data['mag_auto_r'][:] < i+binlen) & (mask_tile))
        
        if len(extent[sliced]) > 30:
            p0 = [5., 2., 0.2]
            amp, mean, stdv = FGtH(extent[sliced], ylims, bb, p0, guess_from_avg=True, show_plot=plot_by_slice)
            shape = stdv/amp
                
            if shape > 0.02: 
                if i > lims[0]:
                    p0 = [10., mean_pb[0], 0.2]
                else :
                    p0 = [10., mean, 0.2]
                amp, mean, stdv = FGtH(extent[sliced], (ylims[0], mean), bb/2, p0, guess_from_avg=False, show_plot=plot_by_slice)

            amp_pb.append(amp)
            mean_pb.append(mean)
            stdv_pb.append(stdv)
            magbins.append(i)
            mean_magerr.append(np.mean(data['rJAVA'][sliced,1]))

    
    magbins = np.array(magbins)
    amp_pb = np.array(amp_pb)
    mean_pb = np.array(mean_pb)
    stdv_pb = np.array(stdv_pb)
    mean_magerr = np.array(mean_magerr)
    
    tot_err = np.empty(len(mean_magerr))
    for j in range(len(mean_magerr)):
        tot_err[j] = np.sqrt(mean_magerr[j]**2. + stdv_pb[j]**2.)
    
    #-- LINEAR FIT --#
    magbin_mask = (magbins <= 19.0)
    a, b = tools.linfit(magbins[magbin_mask], mean_pb[magbin_mask], stdv_pb[magbin_mask])
    
    #-- ERROR FIT --#
    magbin_mask = (magbins <= 20.3)
    par = [1., 0.4, 40.]
    pars, varmatr = curve_fit(magerr_fit, magbins[magbin_mask], tot_err[magbin_mask], p0=par)
    k = pars[0]
    c = pars[1]
    d = pars[2]



    '''
    #-- PLOT --#
    xcs = np.arange(15., 21., 0.05)
    iy = a+b*magbins                                 #linear fit of gaussian centers
    yps = k + (10.**(c*xcs+d)) + (a+b*xcs)           #straight fit to the error
    yps2 = 2*(k + (10.**(c*xcs+d))) + (a+b*xcs)      #translated straight fit
    yps5 = 3*(k + (10.**(c*xcs+d))) + (a+b*xcs)      #translated straight fit
    yps3 = -1.*(k + (10.**(c*xcs+d))) + (a+b*xcs)    #simmetric fit with respect to the linear fit 
    yps4 = mean_magerr/tot_err                       #relative importance of mag error in the total error

    # border = 3*(k + (10.**(c*data['mag_auto_r'][:]+d))) + a + b*data['mag_auto_r'][:]    #decomment to plot the final compact-extended selection
    # extended = ((extent > border) & (mask_tile))                                   #decomment to plot the final compact-extended selection
    # compact = ((extent <= border) & (mask_tile))                                   #decomment to plot the final compact-extended selection
    
    pink = (1.0, 0.6, 0.7)
    pink2 = (1.0, 0.8, 0.9)
    fig = plt.figure(figsize=(12,10))

    s1 = plt.subplot2grid((5, 1), (0, 0), rowspan=3)
    s1.set_ylim( [0., 6.])
    s1.plot(data['mag_auto_r'][mask_tile], extent[mask_tile], 'og', markersize=5, alpha=0.1)     # comment to plot the final compact-extended selection
    # s1.plot(data['mag_auto_r'][extended], extent[extended], 'og', markersize=5, alpha=0.1)     #decomment to plot the final compact-extended selection
    # s1.plot(data['mag_auto_r'][compact], extent[compact], 'oy', markersize=5, alpha=0.1)       #decomment to plot the final compact-extended selection
    s1.plot(magbins, iy, 'k-', label='linear fit of Gaus.centr.val.')
    s1.errorbar(magbins, mean_pb, yerr=stdv_pb, fmt='o', ecolor='r', color='r', markersize=2, linewidth=2, label='Gaussian fit per bin')
    s1.plot(magbins, iy+mean_magerr, 'b-', label='average mag error per bin')
    s1.plot(magbins, iy-mean_magerr, 'b-')
    s1.plot(xcs, yps, 'm-', linewidth=2, label='total-error fit')
    s1.plot(xcs, yps2, color=pink, linewidth=2, label='2 * (total-error fit)')
    s1.plot(xcs, yps5, color=pink2, linewidth=2, label='3 * (total-error fit)')
    s1.plot(xcs, yps3, 'm-', linewidth=2)
    s1.set_title('Tile: ' + str(int(tilenum))+'  morphology')
    s1.set_xlabel('rJAVA  [mags]')
    s1.set_ylabel('mu_max - rJAVA  [auto-mags]')
    s1.legend()

    #---- CLASS_STAR PLOT ----#
    s2 = plt.subplot2grid((5, 1), (3, 0), rowspan=2, sharex=s1)
    s2.plot(data['mag_auto_r'][mask_tile], data['cstar'][mask_tile], 'ob', markersize=2, alpha=0.1)
    s2.set_xlabel('rJAVA  [mags]')
    s2.set_ylabel('CLASS_STAR')

    # #---- PLOT MAG_ERROR "SIGNIFICANCE" within TOTAL ERROR ----#
    # s2 = plt.subplot2grid((5, 1), (4, 0), sharex=s1)
    # s2.plot(magbins, yps4, 'r-', linewidth=2)
    # s2.plot((14., 24.), (0.5, 0.5), 'k--', linewidth=0.75)
    # s2.set_title('mag_err "importance"')
    # s2.set_xlabel('rJAVA  [mags]')
    # s2.set_ylabel('mag/tot err')
    
    plt.tight_layout()
    #plt.savefig(setup['plots']+'interesting_morpho-plot_tile'+str(tilenum)+'_withErrorCurve.png')
    plt.show()
    plt.close()
    '''

    parameters = [a, b , k , c, d]
    return parameters
