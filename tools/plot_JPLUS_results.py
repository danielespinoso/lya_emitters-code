import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as col
import sys
import os

import tools.jplus_filter_system
from tools.converters import magtoflux as MtF
from tools.converters import fluxtomag as FtM

import tools.settings
setup = tools.settings.set_up()

sys.path.append(setup['jplus_code'])
import jplus


def plot_spec(lamb, x, errx, title, unt='mags', limits=[3400., 9000., 25.5, 16.5], idd=0):
    fig = plt.figure(figsize=(12,10))
    matplotlib.rcParams.update({'font.size': 16})
    matplotlib.rcParams.update({'lines.linewidth': 3})
    matplotlib.rcParams.update({'lines.markersize': 6})
    ax = plt.subplot(111)
    kol = (1.,0,1.) #uJAVA
    ax.errorbar(lamb[0], x[0], yerr=errx[0], c=kol, fmt='s', ecolor=kol, alpha=0.7, label='uJAVA')
    kol = (0.7,0,0.7) #J0378
    ax.errorbar(lamb[1], x[1], yerr=errx[1], c=kol, fmt='o', ecolor=kol, alpha=0.7, label='J0378')
    kol = (0.5,0,1.)  #J0395
    ax.errorbar(lamb[2], x[2], yerr=errx[2], c=kol, fmt='o', ecolor=kol, alpha=0.7, label='J0395')
    kol = (0.,0.,1.)  #J0410
    ax.errorbar(lamb[3], x[3], yerr=errx[3], c=kol, fmt='o', ecolor=kol, alpha=0.7, label='J0410')
    kol = (0.,0.8,0.8)  #J0430
    ax.errorbar(lamb[4], x[4], yerr=errx[4], c=kol, fmt='o', ecolor=kol, alpha=0.7, label='J0430')
    kol = (0.,0.6,0.6)  #gJAVA
    ax.errorbar(lamb[5], x[5], yerr=errx[5], c=kol, fmt='s', ecolor=kol, alpha=0.7, label='gJAVA')
    kol = (0.4,0.8,0.)  #J0515
    ax.errorbar(lamb[6], x[6], yerr=errx[6], c=kol, fmt='o', ecolor=kol, alpha=0.7, label='J0515')
    kol = (1.,0.5,0.)  #rJAVA
    ax.errorbar(lamb[7], x[7], yerr=errx[7], c=kol, fmt='s', ecolor=kol, alpha=0.7, label='rJAVA')
    kol = (1.,0.,0.)  #J0660
    ax.errorbar(lamb[8], x[8], yerr=errx[8], c=kol, fmt='o', ecolor=kol, alpha=0.7, label='J0660')
    kol = (0.8,0.,0.)  #iJAVA
    ax.errorbar(lamb[9], x[9], yerr=errx[9], c=kol, fmt='s', ecolor=kol, alpha=0.7, label='iJAVA')
    kol = (0.6,0.,0.)  #J0861
    ax.errorbar(lamb[10], x[10], yerr=errx[10], c=kol, fmt='o', ecolor=kol, alpha=0.7, label='J0861')
    kol = (0.3,0.,0.)  #zJAVA
    ax.errorbar(lamb[11], x[11], yerr=errx[11], c=kol, fmt='s', ecolor=kol, alpha=0.7, label='zJAVA')

    ax.plot(lamb, x, 'k-', alpha=0.3, linewidth=0.8) #solid line to join points

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    # limits=[3400., 9000., 25., 17.]
    ax.axis(limits)
    ax.set_title(title, fontsize=15)
    ax.set_xlabel(r'$\rm{wavelenght}\quad[\AA]$')
    if unt == 'mags':
        ax.set_ylabel(setup['mag_type']+'  [mags]')
    else:
        #ax.set_ylabel(setup['mag_type']+'_flux  [erg/s*cm2]')
        ax.set_ylabel(r'$\rmF_\lambda\quad[erg/s*cm^2]$')
    ax.legend(loc='center left', fontsize = 'small', bbox_to_anchor=(1, 0.5))
    #plt.savefig(setup['plots']+'spectra/spectrum_BPZobj'+str(idd+1)+'_BFtemplate_withANDnoNB.png')
    #plt.savefig(setup['plots']+title+'.png')
    #plt.savefig(setup['plots']+'composite_'+title+'.png')
    #plt.savefig(setup['home']+'Dropbox/lae_ddt_proposal/photospectra/'+title[13:19]+'_photospec.eps', format='eps', dpi=800)
    #plt.savefig(setup['plots']+'meetings/26-04-2016/qso'+str(idd)+'_flux.png')
    #plt.savefig(setup['plots']+setup['mag_type']+'_mags_'+setup['filters'][0]+'_('+setup['data_rels']+'_data)_'+title[:4]+'_candidates_composite.eps', format='eps', dpi=2000)
    #plt.savefig(setup['plots']+setup['mag_type']+'_mags_'+setup['filters'][0]+'_('+setup['data_rels']+'_data)_'+'compact_candidates_median.eps', format='eps', dpi=2000)
    plt.show()
    plt.close()



    

def plot_MEDspectra(data, mask=[], titol=''):
    if len(mask) == 0:
        mask = np.ones(len(data['rJAVA'][:,0]), dtype=bool)
        
    dataset = jplus.tools.select_object(data, mask)
    for i in np.arange(0, 12):
        flt = tools.jpflt(i)
        fl, dfl = MtF(dataset[flt][:,0], band=flt, dmags=dataset[flt][:,1], unit='l')
        key = flt + '_flux'
        dataset[key] = np.array([fl, dfl]).T

    ll = len(dataset['rJAVA'][:,0])

    newspec = {}
    lamb = np.empty(12)
    median = np.empty(12)
    median_err = np.empty(12)
    fref = 1.0e-16  # I choose to put the rJAVA band of all the spectra at 1.e-16
    for i in np.arange(0, 12):
        flt = tools.jpflt(i)
        key = flt + '_flux'
        newspec[flt] = fref*(dataset[key][:,0]/dataset['rJAVA_flux'][:,0])
        lamb[i] = tools.jplus_filter_system.jplus_pivot(flt)
        median[i] = np.median(newspec[flt])
        median_err[i] = np.median(dataset[key][:,1])

    plot_spec(lamb, median, median_err, titol, unt='flux', limits=[3400., 9000., 0., 7.0e-16])
    #return median, median_err


    



# plots the difference between two composite spectra (mask1 - mask2, not the opposite)
def plot_MEDspectra_diff(data, mask1=[], mask2=[], titol=''):
    mask = np.array([mask1, mask2]).T
    lamb = np.empty(12)
    median = np.empty((12,2))
    median_err = np.empty((12,2))
    for qq in range(0,mask.shape[1]) :
        dataset = jplus.tools.select_object(data, mask[:,qq])
        for i in np.arange(0, 12):
            flt = tools.jpflt(i)
            fl, dfl = MtF(dataset[flt][:,0], band=flt, dmags=dataset[flt][:,1], unit='l')
            key = flt + '_flux'
            dataset[key] = np.array([fl, dfl]).T
            
        ll = len(dataset['rJAVA'][:,0])

        newspec = {}
        fref = 1.0e-16  # I choose to put the rJAVA band of all the spectra at 1.e-16 (arbitrary)
        for i in np.arange(0, 12):
            flt = tools.jpflt(i)
            key = flt + '_flux'
            newspec[flt] = fref*(dataset[key][:,0]/dataset['rJAVA_flux'][:,0])
            lamb[i] = tools.jplus_filter_system.jplus_pivot(flt)
            median[i,qq] = np.median(newspec[flt])
            median_err[i,qq] = np.median(dataset[key][:,1])

    diff = median[:,0] - median[:,1]  # ATTENTION!! It's (mask1 - mask2), not the opposite!!
    diff_err = np.sqrt(median_err[:,0]**2. + median_err[:,1]**2.) 
    plot_spec(lamb, diff, diff_err, titol, unt='flux', limits=[3400., 9000., -3.0e-16, 3.0e-16])
    
    #return diff, diff_err

    

    

    
    

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
        lims = [3400., 9000., -0.4e-17, 3.5e-16]
        # LIMITS FOR SDSS QSOs
        # if number == 1:
        #     lims = [3400., 9000., 0.0, 4.0e-16]
        # elif number == 5 or number == 6:
        #     lims = [3400., 9000., 0.0, 4.0e-17]
        # elif number == 11:
        #     lims = [3400., 9000., 0.0, 2.3e-16]
        # elif number == 12:
        #     lims = [3400., 9000., 0.0, 3.0e-16]
        
    #titul = 'ra: '+str(ra)+'  dec: '+str(dec)+'    z: '+str(rsh)
    titul = 'OBJECT:  LAE_extd_1      ra: '+str(ra)+'    dec: '+str(dec)
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


    
def plot_footprint_and_data(jpdata, other=[], plot_jplus=False):
    matplotlib.rcParams.update({'font.size': 14})
    matplotlib.rcParams.update({'lines.markersize': 4})
    fig = plt.figure(figsize=(12,10))
    for i in set(jpdata['tile_id']):
        tilemask = (jpdata['tile_id'] == i)
        minra = min(jpdata['coords'][tilemask,0]); maxra = max(jpdata['coords'][tilemask,0])
        mindec = min(jpdata['coords'][tilemask,1]); maxdec = max(jpdata['coords'][tilemask,1])
        
        #plot jplus
        if plot_jplus == True:
            coordmask = ((jpdata['coords'][:,0] > minra) & (jpdata['coords'][:,0] < maxra) & \
                         (jpdata['coords'][:,1] > mindec) & (jpdata['coords'][:,1] < maxdec))
            plt.plot(jpdata['coords'][coordmask,0], jpdata['coords'][coordmask,1], 'ob', alpha=0.2)

        #plot other data
        if len(other) != 0:
            coordmask = ((other['coords'][:,0] > minra) & (other['coords'][:,0] < maxra) & \
                         (other['coords'][:,1] > mindec) & (other['coords'][:,1] < maxdec))
            plt.plot(other['coords'][coordmask,0], other['coords'][coordmask,1], 'or', alpha=0.2)

        #plot jplus footprint (one square per tile)
        plt.plot( (minra, maxra), (mindec, mindec), 'k-', alpha=0.8 )
        plt.plot( (maxra, maxra), (mindec, maxdec), 'k-', alpha=0.8 )
        plt.plot( (maxra, minra), (maxdec, maxdec), 'k-', alpha=0.8 )
        plt.plot( (minra, minra), (maxdec, mindec), 'k-', alpha=0.8 )

    
    plt.title(setup['data_rels']+' footprint')
    plt.xlabel('ra [deg]')
    plt.ylabel('dec [deg]')
    plt.legend()
    #plt.show()
    #plt.close()

    

def composite_spec(data, mask=[], titol=''):
    if len(mask) == 0:
        mask = np.ones(len(data['rJAVA'][:,0]), dtype=bool)

    g = tools.singleMtoF(data['gJAVA'][0,0], band='gJAVA')
    r = tools.singleMtoF(data['rJAVA'][0,0], band='rJAVA')
    fref = (g+r)/2. # I choose to normalize all the spectra to the mean of rJAVA-gJAVA bands of the first spectra

    mask = ((mask) & (data['J0395'][:,0] > 18.8))
    dataset = jplus.tools.select_object(data, mask)

    # construct a matrix with fluxes and flux-errors
    for i in np.arange(0, 12):
        flt = tools.jpflt(i)
        fl, dfl = MtF(dataset[flt][:,0], band=flt, dmags=dataset[flt][:,1], unit='l')
        key = flt + '_flux'
        dataset[key] = np.array([fl, dfl]).T

    ll = len(dataset['rJAVA'][:,0])

    newspec = {}
    ww = {}
    lamb = np.empty(12)
    compo = np.empty(12)
    compo_err = np.empty(12)
    # fref = 1.0e-16  # I choose to put the rJAVA band of all the spectra at 1.e-16
    for i in np.arange(0, 12):
        flt = tools.jpflt(i)
        key = flt + '_flux'
        # newspec[flt] = fref*(dataset[key][:,0]/dataset['rJAVA_flux'][:,0])
        DD = (dataset['rJAVA_flux'][:,0]+dataset['gJAVA_flux'][:,0])/2.
        newspec[flt] = fref*(dataset[key][:,0]/DD)
        ww[flt] = (1./dataset[key][:,1])
        mask = (ww[flt] == -1.*np.inf)
        ww[flt][mask] = 0.
        lamb[i] = tools.jplus_filter_system.jplus_pivot(flt)
        compo[i] = np.average(newspec[flt], weights=ww[flt])
        compo_err[i] = np.std(newspec[flt])/np.sqrt(len(newspec[flt]))#np.median(dataset[key][:,1])
        # if compo_err[i] > 4.e-18:
        #     compo_err[i] = 1.5e-18
        #compo_err[i] = np.sqrt(sum(dataset[key][:,1]**2.))/np.sqrt(len(dataset[key][:,1]))#np.median(dataset[key][:,1])

    tools.plot_spec(lamb, compo, compo_err, titol, unt='flux', limits=[3400., 9000., 0., 1.e-16])






def plot_SDSSandJPLUS_spectra(sdss_file, first_data, first_objid, mask=[], units='mags', zsdss=0, zfromSDSS=False, number=0):
    # SDSS SPECTRUM
    # wave, flu = np.genfromtxt(sdss_file, unpack=True, usecols=(0,1), delimiter=',', skip_header=1)

    import pyfits as pf
    data, header = pf.getdata(sdss_file, 0, header=True)
    flu = np.empty(len(data))
    wave = np.empty(len(data))
    for i,j in zip(range(len(data)), data):
        flu[i] = j[0]
        wave[i] = 10**j[1]
    
    flu = flu*1.3e-17   # normalization to match jplus data. SDSS spectra are in units of 10^-17 erg/s/cm2/A
    
    f=plt.figure(figsize=(12,10))
    plt.plot(wave, flu, 'k-', linewidth=0.5, alpha=0.3, label='SDSS')

    # JPLUS SPECTRUM
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
        lims = [3400., 9000., -5.0e-17, 3.3e-16]

    titul = 'object: SDSS QSO     ra: '+str(ra)+'   dec: '+str(dec)+'     z: '+str(rsh)
        
    matplotlib.rcParams.update({'font.size': 18})
    matplotlib.rcParams.update({'lines.linewidth': 5})
    matplotlib.rcParams.update({'lines.markersize': 10})
    
    kol = (1.,0,1.) #uJAVA
    plt.errorbar(lamb[0], values[0], yerr=errors[0], c=kol, fmt='s', ecolor=kol, alpha=0.95, label='uJAVA')
    kol = (0.7,0,0.7) #J0378
    plt.errorbar(lamb[1], values[1], yerr=errors[1], c=kol, fmt='s', mfc='white', ecolor=kol, alpha=0.95, label='J0378')
    kol = (0.5,0,1.)  #J0395
    plt.errorbar(lamb[2], values[2], yerr=errors[2], c=kol, fmt='s', mfc='white', ecolor=kol, alpha=0.95, label='J0395')
    kol = (0.,0.,1.)  #J0410
    plt.errorbar(lamb[3], values[3], yerr=errors[3], c=kol, fmt='s', mfc='white', ecolor=kol, alpha=0.95, label='J0410')
    kol = (0.,0.8,0.8)  #J0430
    plt.errorbar(lamb[4], values[4], yerr=errors[4], c=kol, fmt='s', mfc='white', ecolor=kol, alpha=0.95, label='J0430')
    kol = (0.,0.6,0.6)  #gJAVA
    plt.errorbar(lamb[5], values[5], yerr=errors[5], c=kol, fmt='s', ecolor=kol, alpha=0.95, label='gJAVA')
    kol = (0.4,0.8,0.)  #J0515
    plt.errorbar(lamb[6], values[6], yerr=errors[6], c=kol, fmt='s', mfc='white', ecolor=kol, alpha=0.95, label='J0515')
    kol = (1.,0.5,0.)  #rJAVA
    plt.errorbar(lamb[7], values[7], yerr=errors[7], c=kol, fmt='s', ecolor=kol, alpha=0.95, label='rJAVA')
    kol = (1.,0.,0.)  #J0660
    plt.errorbar(lamb[8], values[8], yerr=errors[8], c=kol, fmt='s', mfc='white', ecolor=kol, alpha=0.95, label='J0660')
    kol = (0.8,0.,0.)  #iJAVA
    plt.errorbar(lamb[9], values[9], yerr=errors[9], c=kol, fmt='s', ecolor=kol, alpha=0.95, label='iJAVA')
    kol = (0.6,0.,0.)  #J0861
    plt.errorbar(lamb[10], values[10], yerr=errors[10], c=kol, fmt='s', mfc='white', ecolor=kol, alpha=0.95, label='J0861')
    kol = (0.3,0.,0.)  #zJAVA
    plt.errorbar(lamb[11], values[11], yerr=errors[11], c=kol, fmt='s', ecolor=kol, alpha=0.95, label='zJAVA')

    plt.axis(lims)
    plt.title(titul, fontsize=18)
    plt.xlabel(r'$\rm{wavelenght}\quad[\AA]$', fontsize=18)
    if units == 'mags':
        plt.ylabel(setup['mag_type']+'  [mags]', fontsize=18)
    else:
        plt.ylabel(r'$\rmF_\lambda\quad[\,erg\,/\,(\AA\,s\,cm^2)\,]$', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.legend(loc='upper right', fontsize = 'small')
    plt.savefig(setup['plots']+'JPLUS_EDRobj-8461-19884_and_SDSS_spec-4447-55542-0346.pdf')
    #plt.show()
    plt.close()
