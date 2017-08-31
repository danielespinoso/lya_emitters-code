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
def evaluate_morphology(data, extent, mask_tile, tilenum, plot_by_slice = False, show_plt=False):
    from tools.morphology_maglimits_dictionary import maglimits
    lims = maglimits(tilenum)
    binlen = 0.25
    binnum = 50
    fr_lim = [14.5, 15.0]                         # first-reference by-eye limits
    if (tilenum > 10015): fr_lim = [16.5, 17.0]   # first-reference by-eye limits   

    amp_pb = []       #amplitude of the first gaussian fit
    mean_pb = []      #first-gaussian central-value at each mag_bin
    stdv_pb = []      #standard deviation of the first gaussian distribution
    magbins = []      #mag vector to plot first gaussian fits
    mean_magerr = []  #mean error on magnitudes as a function of magnitudes
    
    sliced = ( (data['mag_auto_r'] > fr_lim[0]) & (data['mag_auto_r'] < fr_lim[1]) & (mask_tile) )
    if tilenum == 2521:
        sliced = ( (data['mag_auto_r'] > 17.) & (data['mag_auto_r'] < 17.5) & (mask_tile) )
    mdn = np.median(extent[sliced])   # now we have a tile-dependent reference for the "safe_box"
    hlims = [mdn-0.55, mdn+0.55]
    safe_box = ( (data['mag_auto_r'] > lims[0]) & (data['mag_auto_r'] < lims[1]) & \
                 (extent > hlims[0]) & (extent < hlims[1]) & (mask_tile) )

    # # TEST PLOT - BEFORE ALL OPERATIONS
    # plt.plot(data['mag_auto_r'][mask_tile], extent[mask_tile], 'og', alpha=0.15, markersize=2)
    # plt.plot(data['mag_auto_r'][safe_box], extent[safe_box], 'ob', alpha=0.6, markersize=2)
    # plt.show()
    # plt.close()
    # sys.exit()
    

    for i in np.arange(lims[0], lims[1], binlen):
        sliced = ((data['mag_auto_r'][:] > i) & (data['mag_auto_r'][:] < i+binlen) & (mask_tile))
        
        if len(extent[sliced]) > 15:
            p0 = [5., 2., 0.2]
            amp, mean, stdv = FGtH(extent[sliced], hlims, binnum, p0, guess_from_avg=True, show_plot=plot_by_slice)
            
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
    a, b = tools.linfit(magbins, mean_pb, stdv_pb)
    
    #-- ERROR FIT --#
    par = [0.1, 0.4, -10.]
    if (tilenum == 10053) or (tilenum == 9859) or (tilenum == 11365) or (tilenum == 11367) or (tilenum == 11371) or (tilenum == 10978):
        pars = [0.1, 0.4, -10.]
    else:
        pars, varmatr = curve_fit(magerr_fit, magbins, tot_err, p0=par)
    k = pars[0]
    c = pars[1]
    d = pars[2]

    

    #------- THE BIG PLOT -------#
    if show_plt == True:
        import matplotlib
        matplotlib.rcParams.update({'font.size': 16})
        matplotlib.rcParams.update({'font.weight': 'bold'})
        matplotlib.rcParams.update({'lines.linewidth': 2})
        matplotlib.rcParams.update({'lines.markersize': 4})
        pink = (1.0, 0.6, 0.7)
        pink2 = (1.0, 0.8, 0.9)
        
        # # LINES TO PLOT
        xcs = np.arange(13., 23., 0.05)
        iy = a+b*magbins                                 #linear fit of gaussian centers
        yps = k + (10.**(c*xcs+d)) + (a+b*xcs)           #straight fit to the error
        yps2 = 2*(k + (10.**(c*xcs+d))) + (a+b*xcs)      #translated straight fit
        yps5 = 3*(k + (10.**(c*xcs+d))) + (a+b*xcs)      #translated straight fit
        yps3 = -1.*(k + (10.**(c*xcs+d))) + (a+b*xcs)    #simmetric fit with respect to the linear fit 
        yps4 = mean_magerr/tot_err                       #relative importance of mag error in the total error
        # border = 3*(k + (10.**(c*data['mag_auto_r'][:]+d))) + a + b*data['mag_auto_r'][:]    #decomment to plot the final compact-extended selection
        # extended = ((extent > border) & (mask_tile))                                         #decomment to plot the final compact-extended selection
        # compact = ((extent <= border) & (mask_tile))                                         #decomment to plot the final compact-extended selection
        
        # # ACTUAL PLOT
        fig = plt.figure(figsize=(12,10))

        # GUILLAUME MORPHOLOGY SPACE
        s1 = plt.subplot2grid((5, 1), (0, 0), rowspan=3)
        s1.set_xlim( [13., 24.5])
        s1.set_ylim( [0., 6.])
        s1.plot(data['mag_auto_r'][mask_tile], extent[mask_tile], 'og', alpha=0.08)     # comment to plot the final compact-extended selection
        # s1.plot(data['mag_auto_r'][extended], extent[extended], 'og', alpha=0.1)     #decomment to plot the final compact-extended selection
        # s1.plot(data['mag_auto_r'][compact], extent[compact], 'oy', alpha=0.1)       #decomment to plot the final compact-extended selection
        s1.errorbar(magbins, mean_pb, yerr=stdv_pb, fmt='o', ecolor='r', color='r', linewidth=4, label='Fit per bin')
        s1.plot(magbins, iy, 'k-', label='linear fit')
        s1.plot(magbins, iy+mean_magerr, 'c-', label=r'$\sigma_{NB}$'+' per bin')
        s1.plot(magbins, iy-mean_magerr, 'c-')
        #s1.plot(xcs, yps, color=pink2, label='Error fit')          # 1-sigma total-error fit (up)   -   light pink
        #s1.plot(xcs, yps3, color=pink2)                            # 1-sigma total-error fit (down) -   light pink
        #s1.plot(xcs, yps2, color=pink, label='2 * (Error fit)')    # 2-sigma total-error fit        -   dark pink
        s1.plot(xcs, yps5, 'm-', label='Threshold')                 # 3-sigma total-error fit        -   magenta
        s1.set_title('Tile: ' + str(int(tilenum))+'  morphology', fontweight='bold')
        s1.set_ylabel('mu_max - rJAVA  [auto-mags]', fontweight='bold')
        s1.text(14., 5., '"Extended region"', color='m')
        s1.text(18., 0.8, '"Compact region"', color='m')
        s1.legend(fontsize=14)
        
        # CLASS_STAR - MAGNITUDE SPACE 
        s2 = plt.subplot2grid((5, 1), (3, 0), rowspan=2, sharex=s1)
        s2.set_xlim( [13., 24.5])
        s2.plot(data['mag_auto_r'][mask_tile], data['cstar'][mask_tile], 'ob', markersize=2, alpha=0.1)
        s2.plot((13., 24.5), (0.9, 0.9), 'r--')
        s2.text(22., 0.95, '"Extended region"', color='r', fontsize=10)
        s2.text(22., 0.77, '"Compact region"', color='r', fontsize=10)
        s2.text(14., 0.8, 'CUT = 0.9', color='r', fontsize=10)
        s2.set_xlabel('rJAVA  [mags]', fontweight='bold')
        s2.set_ylabel('CLASS_STAR', fontweight='bold')
        
        # # MAG_ERROR "SIGNIFICANCE" within TOTAL ERROR
        # s2 = plt.subplot2grid((5, 1), (4, 0), sharex=s1)
        # s2.plot(magbins, yps4, 'r-', linewidth=2)
        # s2.plot((14., 24.), (0.5, 0.5), 'k--', linewidth=0.75)
        # s2.set_title('mag_err "importance"')
        # s2.set_xlabel('rJAVA  [mags]')
        # s2.set_ylabel('mag/tot err')
        
        plt.tight_layout()
        fold = 'CLASS_STAR_comparison_and_PROPER_morpho_criterion_(gauss_fit+error_fit+3_sigma_line)/'
        #plt.savefig(setup['plots']+'T3_morphoplots/'+'morpho-plot_tile'+str(tilenum)+'_cstar_errThreshold.png')
        #plt.savefig(setup['plots']+'morphology/mu_max/'+fold+'morpho-plot_tile'+str(tilenum)+'_cstar_errThreshold.eps', format='eps', dpi=1000)
        #plt.savefig(setup['plots']+'morphology/mu_max/'+fold+'morpho-plot_tile'+str(tilenum)+'_cstar_errThreshold.pdf')
        plt.show()
        plt.close()
    

    parameters = [a, b , k , c, d]
    return parameters



















'''
###################################
#                                 #
#     OLD MORPHOLOGY ANALYSIS     #
#                                 #
###################################

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
                    #p0 = [10., mean_pb[0], 0.2]
                    p0 = [10., mean, 0.2]
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
    

    parameters = [a, b , k , c, d]
    return parameters
'''
