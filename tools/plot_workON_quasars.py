import numpy as np
import sys

import matplotlib.pyplot as plt

import tools.settings
setup = tools.settings.set_up()

def plot_quasar_ColorColor(jplus, quasars, mask_quasars, index):
    #-- COLOR-COLOR DIAGRAMS (QSO_vs_JPLUS) --# 
    xcol = ['uJAVA', 'iJAVA']
    ycol = ['gJAVA', 'zJAVA']
    tagl = 22.5

    cmask = ((jplus[xcol[0]][:,0] > 0) & (jplus[xcol[0]][:,0] < tagl) &\
             (jplus[xcol[1]][:,0] > 0) & (jplus[xcol[1]][:,0] < tagl) &\
             (jplus[ycol[0]][:,0] > 0) & (jplus[ycol[0]][:,0] < tagl) &\
             (jplus[ycol[1]][:,0] > 0) & (jplus[ycol[1]][:,0] < tagl) )
    qsocmask = ((mask_quasars) & (cmask))
    jplcmask = ((~mask_quasars) & (cmask))
    zmask = (quasars['redshift'][index[qsocmask]] > 2.0)
    jpl_col_x = jplus[xcol[0]][jplcmask,0] - jplus[xcol[1]][jplcmask,0]
    jpl_col_y = jplus[ycol[0]][jplcmask,0] - jplus[ycol[1]][jplcmask,0]
    qso_col_x = jplus[xcol[0]][qsocmask,0] - jplus[xcol[1]][qsocmask,0]
    qso_col_y = jplus[ycol[0]][qsocmask,0] - jplus[ycol[1]][qsocmask,0]
    
    plt.plot(jpl_col_x, jpl_col_y, 'ob', markersize=4, alpha=0.4, label='not QSO')
    plt.plot(qso_col_x[~zmask], qso_col_y[~zmask], 'or', markersize=4, alpha=0.4, label='z<2 QSO')
    plt.plot(qso_col_x[zmask], qso_col_y[zmask], 'og', markersize=4, alpha=0.4, label=' z>2 QSO')
    plt.xlabel(xcol[0]+' - '+xcol[1])
    plt.ylabel(ycol[0]+' - '+ycol[1])
    plt.legend()
    #plt.savefig(setup['plots']+'UG-GR_color-color_QSOvsCANDIDATES.png')
    plt.show()
    plt.close()



    
def plot_quasar_zDistrib(quasars, mask_quasars, index):
    #-- REDSHIFT DISTRIBUTION of SDSS QSOs --#
    plt.hist(quasars['redshift'][index[mask_quasars]], range=(0.0, 2.5),\
             bins=25, edgecolor='k', alpha=0.6)
    plt.xlabel('z')
    plt.title('x-matched SDSS QSOs z distribution')
    plt.show()
    #plt.savefig(setup['plots']+'sdssQSO_contaminants_zdistrib.png')
    plt.close()




def plot_quasar_PhotoSpec(jplus, quasars, mask_quasars, index):
    from tools.plot_JPLUS_results import plot_JPLUSphotoSpectra as PhotoSpeck
    #-- SDSS QSO Photo-spectra (According to jplus photometry) --#
    cc = 0
    zz = quasars['redshift'][index[mask_quasars]]     # redshift from SDSS
    for i in np.arange(len(jplus['coords'][:,0])):
        if mask_quasars[i] == True:
            PhotoSpeck(jplus, i, units='flux', zsdss=zz[cc], zfromSDSS=True, number=cc+1)
            cc += 1



def plot_quasar_lumifunc(jplus, quasars, mask_quasars, index, n_tiles):
    from scipy.interpolate import interp1d
    from tools.threeFM import jplus3FM
    
    #-- Lya LUMINOSITY DISTRIBUTION of SDSS QSOs --#
    fluxlya, fcont, flux_nb, fcont_error = jplus3FM(jplus)
    fluxlya = fluxlya[mask_quasars]    # selects only SDSS QSOs
    comov_table = setup['tools'] + 'comoving_distance.txt'
    z, eta = np.genfromtxt(comov_table, skip_header=1, unpack=True)
    line = interp1d(z, eta, kind='nearest', bounds_error=None, fill_value='extrapolate')
    zz = quasars['redshift'][index[mask_quasars]]     # redshift from SDSS
    dl = line(zz)*(1+zz)      #now in Mpc, redshift from SDSS
    dl = dl*(1.e6)*3.0856e18              #now in cm 
    lumilya = 4.*np.pi*fluxlya*(dl**2.)   #now in erg/s
    loglum = np.log10(lumilya)
    # plt.hist(loglum, range=(43, 45.5), color='r', bins=9, edgecolor='r', alpha=0.6)
    # plt.xlabel('Log(Lya)')
    # plt.title('Lya luminosity distribution SDSS QSOs z ')
    # plt.show()
    # plt.close()

    #-- TENTATIVE LUMINOSITY FUNCTION OF SDSS QSO in JPLUS candidates --#
    jplus_area = (1.4**2.)*n_tiles #now in sq.degrees
    sphere_area = 4*np.pi/((np.pi/180.)**2.) #now in sq.degrees
    fsky = jplus_area/sphere_area
    dc_min = line(min(zz))   #now in Mpc
    dc_max = line(max(zz))   #now in Mpc
    dVc = (4.*np.pi*(dc_max**3. - dc_min**3.)/3.)*fsky
    dV = dVc
    dLogL = 0.35
    from tools.histo import histogram
    Lmin = min(loglum) ; Lmax = max(loglum)
    centers, histo = histogram(loglum, Lmin, Lmax, dLogL)
    phi = histo*(1./(dV*dLogL))
    errs = (1./np.sqrt(histo))*(1./(dV*dLogL))
    plt.errorbar(centers, phi, yerr=errs, fmt='or', ecolor='b', markersize=3, label='JPLUS data')
    plt.xlabel('Log L  [erg/s]', fontsize=12)
    plt.ylabel(r'$\Phi\quad  1\,/\,[\rm{Mpc}^3\ \rm{dLog L}]$', fontsize=12)
    plt.yscale('log')
    plt.title('Luminosity function  -  SDSS QSOs in lya candidates')
    #plt.legend(loc=4)
    plt.show()
    plt.close()
    
