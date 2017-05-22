import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import sys
import os

import tools.jplus_filter_system
import tools.converters

import tools.settings
setup = tools.settings.set_up()

sys.path.append(setup['jplus_code'])
import jplus


def plot_spec(lamb, x, errx, title, unt='mags', limits=[3400., 9000., 25.5, 16.5], idd=0):
    fig = plt.figure()
    ax = plt.subplot(111)
    kol = (1.,0,1.) #uJAVA
    ax.errorbar(lamb[0], x[0], yerr=errx[0], c=kol, fmt='s', ecolor=kol, markersize=4, alpha=0.7, label='uJAVA')
    kol = (0.7,0,0.7) #J0378
    ax.errorbar(lamb[1], x[1], yerr=errx[1], c=kol, fmt='o', ecolor=kol, markersize=4, alpha=0.7, label='J0378')
    kol = (0.5,0,1.)  #J0395
    ax.errorbar(lamb[2], x[2], yerr=errx[2], c=kol, fmt='o', ecolor=kol, markersize=4, alpha=0.7, label='J0395')
    kol = (0.,0.,1.)  #J0410
    ax.errorbar(lamb[3], x[3], yerr=errx[3], c=kol, fmt='o', ecolor=kol, markersize=4, alpha=0.7, label='J0410')
    kol = (0.,0.8,0.8)  #J0430
    ax.errorbar(lamb[4], x[4], yerr=errx[4], c=kol, fmt='o', ecolor=kol, markersize=4, alpha=0.7, label='J0430')
    kol = (0.,0.6,0.6)  #gJAVA
    ax.errorbar(lamb[5], x[5], yerr=errx[5], c=kol, fmt='s', ecolor=kol, markersize=4, alpha=0.7, label='gJAVA')
    kol = (0.4,0.8,0.)  #J0515
    ax.errorbar(lamb[6], x[6], yerr=errx[6], c=kol, fmt='o', ecolor=kol, markersize=4, alpha=0.7, label='J0515')
    kol = (1.,0.5,0.)  #rJAVA
    ax.errorbar(lamb[7], x[7], yerr=errx[7], c=kol, fmt='s', ecolor=kol, markersize=4, alpha=0.7, label='rJAVA')
    kol = (1.,0.,0.)  #J0660
    ax.errorbar(lamb[8], x[8], yerr=errx[8], c=kol, fmt='o', ecolor=kol, markersize=4, alpha=0.7, label='J0660')
    kol = (0.8,0.,0.)  #iJAVA
    ax.errorbar(lamb[9], x[9], yerr=errx[9], c=kol, fmt='s', ecolor=kol, markersize=4, alpha=0.7, label='iJAVA')
    kol = (0.6,0.,0.)  #J0861
    ax.errorbar(lamb[10], x[10], yerr=errx[10], c=kol, fmt='o', ecolor=kol, markersize=4, alpha=0.7, label='J0861')
    kol = (0.3,0.,0.)  #zJAVA
    ax.errorbar(lamb[11], x[11], yerr=errx[11], c=kol, fmt='s', ecolor=kol, markersize=4, alpha=0.7, label='zJAVA')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.axis(limits)
    ax.set_title(title)
    ax.set_xlabel('wavelenght  [A]')
    if unt == 'mags':
        ax.set_ylabel('magnitude  [mags]')
    else:
        ax.set_ylabel('flux  [erg/s*cm2]')
    ax.legend(loc='center left', fontsize = 'small', bbox_to_anchor=(1, 0.5))
    #plt.savefig(setup['plots']+'spectra/spectrum_BPZobj'+str(idd+1)+'_BFtemplate_withANDnoNB.png')
    #plt.savefig(setup['plots']+'meetings/26-04-2016/qso'+str(idd)+'_flux.png')
    plt.show()
    plt.close()



    

def plot_MEDspectra(data, mask):
    lamb = np.empty(12)
    mags = np.empty(12)
    errmags = np.empty(12)
    fluxs = np.empty(12)
    errfluxs = np.empty(12)
    X = [],[],[],[],[],[],[],[],[],[],[],[]
    for i,j in zip(np.arange(0,12), jplus.datasets.jplus_filter_names(only_bb=False) ):
        lamb[i] = tools.jplus_filter_system.jplus_pivot(j)
        aa = data[j][mask,0]
        mm = ( aa > 0.)
        ff = tools.converters.magtoflux(aa[mm], band=j)
        X[i].extend(ff)
        fluxs[i] = np.median(X[i])
        errfluxs[i] = np.std(X[i])
        mags[i], errmags[i] =  tools.converters.singleFtoM(fluxs[i], band=j, dflux=errfluxs[i])

    titul = 'Tentative composite spectra'
    plot_spec(lamb, mags, errmags, titul, unt='mags', limits=[3400., 9000., 23., 19.])
    plot_spec(lamb, fluxs, errfluxs, titul, unt='flux', limits=[3400., 9000., -2.e-17, 1.5e-16])



    

def plot_JPLUSphotoSpectra(first_data, first_objid, mask=[], units='mags', zsdss=0, zfromSDSS=False, number=0):
    if len(mask) == 0:
        mask = np.ones(len(first_data['rJAVA'][:,0]), dtype=bool)
        objid = first_objid
    else:
        objid = int(len(first_data['rJAVA'][:,0])*first_objid/len(first_data['rJAVA'][:,0]))
        
    data = jplus.tools.select_object(first_data, mask)
    
    lamb = np.empty(12)
    mag = np.empty(12); flux = np.empty(12)
    errmag = np.empty(12); errflux = np.empty(12)
    for i, ifilter in zip(np.arange(0,12), jplus.datasets.jplus_filter_names(only_bb=False)):
        lamb[i] = tools.jplus_filter_system.jplus_pivot(ifilter)
        mag[i] = data[ifilter][objid,0]
        errmag[i] = data[ifilter][objid,1]
        flux[i], errflux[i] = tools.singleMtoF(mag[i], band=ifilter, dmags=errmag[i], unit='l')
        

    ra = (int(data['coords'][objid,0]*10000.))/10000.
    dec = (int(data['coords'][objid,1]*10000.))/10000.
    if zfromSDSS == True:
        rsh = (int(zsdss*1000.))/1000.
    else:
        rsh = (int(data['redshift'][objid]*1000.))/1000.

    if units == 'mags' :
        values = mag
        errors = errmag
        lims=[3400., 9000., 25.5, 16.5]
    elif units == 'flux' :
        values = flux
        errors = errflux
        lims = [3400., 9000., -0.4e-17, 0.5e-16]
        # LIMITS FOR SDSS QSOs
        # if number == 1:
        #     lims = [3400., 9000., 0.0, 4.0e-16]
        # elif number == 5 or number == 6:
        #     lims = [3400., 9000., 0.0, 4.0e-17]
        # elif number == 11:
        #     lims = [3400., 9000., 0.0, 2.3e-16]
        # elif number == 12:
        #     lims = [3400., 9000., 0.0, 3.0e-16]
        
    titul = 'ra: '+str(ra)+'  dec: '+str(dec)+'    z: '+str(rsh)
    plot_spec(lamb, values, errors, titul, unt=units, limits=lims, idd=number)


    #-------------- PLOT ALL BPZ TEMPLATES --------------#
    # fig = plt.figure()
    # ax = plt.subplot(111)
    # for i,j in zip(tmpl_list, np.arange(len(tmpl_list))):
    #     rf_lamb, tmplt = np.genfromtxt(home+'BPZ/PSFsimulations/pros/bpz/SED/'+i, unpack=True)
    #     obs_lamb = rf_lamb*(1.+ 2.24)  all templates at z=2.24
    #     shift_tmplt = tmplt*(6.**j)
    #     mask = ((obs_lamb > 2800.) & (obs_lamb < 9000.))
    #     tem_kol = ((11-j)/11., 0., j/11.)
    #     ax.plot(obs_lamb[mask], shift_tmplt[mask], c=tem_kol, linewidth=3, alpha=0.8, label=i)
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    # ax.set_title('All BPZ templates at z=2.24')
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # ax.set_xlabel('wavelenght  [A]')
    # ax.set_ylabel(r'$f_{\lambda}\quad\rm{[erg/s cm}^2\rm{A}]$')
    # ax.set_yscale('log')
    # plt.show()




    
def plot_lumiFunc(fsky, comovDist_interpolation, properVol=False):
    lum = np.genfromtxt(setup['final_list'], skip_header=1, usecols=(5), unpack=True)
    Loglum = np.log10(lum) ;  notinf = (Loglum > 0) ;  Loglum = Loglum[notinf]

    if properVol == True:
        dp_min = comovDist_interpolation(setup['zmin'])/(1+setup['zmin'])   #now in Mpc
        dp_max = comovDist_interpolation(setup['zmax'])/(1+setup['zmax'])   #now in Mpc
        dVp = abs((4.*np.pi*(dp_max**3. - dp_min**3.)/3.))*fsky
        dV = dVp
    else:
        dc_min = comovDist_interpolation(setup['zmin'])   #now in Mpc
        dc_max = comovDist_interpolation(setup['zmax'])   #now in Mpc
        dVc = (4.*np.pi*(dc_max**3. - dc_min**3.)/3.)*fsky
        dV = dVc

    if setup['galexmask'] == True:
        dLogL = 0.12
    else:
        dLogL = 0.15

    from tools.histo import histogram
    Lmin = min(Loglum) ; Lmax = max(Loglum)
    centers, histo = histogram(Loglum, Lmin, Lmax, dLogL)
    phi = histo*(1./(dV*dLogL))
    errs = (1./np.sqrt(histo))*(1./(dV*dLogL))

    kon_list = setup['home'] + 'works/lya_emitters/datasets/konno_points.txt'
    kon_Loglum, kon_phi, kon_errs = np.genfromtxt(kon_list,skip_header=1,unpack=True,usecols=(0,1,2))
    
    plt.errorbar(centers, phi, yerr=errs, fmt='or', ecolor='b', markersize=3, label='JPLUS data')
    plt.errorbar(kon_Loglum, kon_phi, yerr=kon_errs, fmt='og', ecolor='m', markersize=3, label='Konno et al.')
    plt.xlabel('Log L  [erg/s]', fontsize=12)
    plt.ylabel(r'$\Phi\quad  1\,/\,[\rm{Mpc}^3\ \rm{dLog L}]$', fontsize=12)
    plt.yscale('log')
    plt.title('Luminosity function  -  JPLUS lya candidates')
    #plt.savefig(setup['plots'] + 'LumiFunc_z2.24_'+str(dLogL)+'bin_better.png')
    plt.legend()
    plt.show()
    plt.close()


def morpHisto(jpl, sdss, qso, gaia):
    violet = (0.6, 0.0, 0.6)
    orange = (1.0, 0.5, 0.0)
    yellow = (1.0, 1.0, 0.0)
    
    from tools.histo import histogram
    jpl_cent, jpl_hist = histogram(jpl, 0, 10, 0.25)
    sdss_cent, sdss_hist = histogram(sdss, 0, 10, 0.25)
    qso_cent, qso_hist = histogram(qso, 0, 10, 0.25)
    gaia_cent, gaia_hist = histogram(gaia, 0, 10, 0.25)
    
    jpl_hist = jpl_hist/max(jpl_hist)
    sdss_hist = sdss_hist/max(sdss_hist)
    qso_hist = qso_hist/max(qso_hist)
    gaia_hist = gaia_hist/max(gaia_hist)

    plt.plot(jpl_cent, jpl_hist, 'b', alpha=0.4, label='selected')
    plt.plot(sdss_cent, sdss_hist, color=violet, markersize=6, alpha=0.4, label='galaxies')
    plt.plot(qso_cent, qso_hist, color=yellow, markersize=6, alpha=0.8, label='qso')
    plt.plot(gaia_cent, gaia_hist, color='r', markersize=6, alpha=0.4, label='gaia')
    plt.xlabel(' FWHM / PDF ')
    plt.title('FWHM and PSF both in rJAVA')
    plt.legend()
    #plt.savefig(setup['plots']+'morpho_hist.png')
    plt.show()
    plt.close()
