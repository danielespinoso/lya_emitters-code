import sys
import time
import deepdish as dd
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d

import tools

setup = tools.set_up()

sys.path.append(setup['jplus_code'])
import jplus
from jplus.tools import crossmatch_angular as xmatch



#######################################################
#
#                WORK ON JPLUS CANDIDATES
#
#######################################################


#--------------- LOADING DATASETS ---------------#
jpl = dd.io.load(setup['jplus_candidates'])       # loads the list of candidates
jpl['redshift'] = np.zeros(len(jpl['coords'][:,0]))
#alljpl = dd.io.load(setup['jplus_input'])        # loads the initial jplus catalogue
    
if setup['load_sdssGal'] == True:
    sdss_gal = dd.io.load(setup['sdss_gal'])      #load sdss galaxies catalogue

if setup['load_sdssQSO'] == True:
    qso = dd.io.load(setup['sdss_qso'])           #load sdss quasars catalogue
    qso_zmask = ((qso['zspec'] > setup['zmin']) & (qso['zspec'] < setup['zmax']))
    qso_z = jplus.tools.select_object(qso, qso_zmask)
    
if setup['load_sdssSTAR'] == True:
    sdss_str = dd.io.load(setup['sdss_str'])      #load sdss stars catalogue
    
if setup['load_gaia'] == True:
    gaia = dd.io.load(setup['gaia'])              #load gaia stars catalogue

if setup['load_mock_cd'] == True:
    mock = dd.io.load(setup['mock_candidates'])   #load mocks catalogue
    
if setup['galexmask'] == True:
    galex = dd.io.load(setup['galex_in'])         #load galex catalogue
    
if setup['load_rafa'] == True:
    rafa = dd.io.load(setup['rafa'])              #load rafa catalog of HII regions in jplus
    
if setup['load_lqac'] == True:
    lqac = dd.io.load(setup['lqac'])              #load Large Qasar Astrometric Catalogue
    
print '\n',time.strftime("%d/%m/%Y")
print '---------- GENERAL INFO ----------'
print 'Filter: ', setup['filters'][0]
print 'Mag type: '+setup['mag_type']+' mags'
print '\nList of loaded datasets:'
if setup['load_sdssGal'] == True: print 'sdss galaxies'
if setup['load_sdssQSO'] == True: print 'sdss quasars'
if setup['load_sdssSTAR'] == True: print 'sdss stars'
if setup['load_gaia'] == True: print 'overlapping gaia'
if setup['load_mock_cd'] == True: print 'jplus mocks'
if setup['galexmask'] == True: print 'galex sources'
if setup['load_rafa'] == True: print 'rafa HII sources'
if setup['load_lqac'] == True: print 'lqac catalogue'

print '\nInitial candidates: ', len(jpl['coords'][:,0])





#-------------- BB SIGNAL-to-NOISE SELECTION  --------------#
# In general, this selection MUST be done tile-by-tile (JPLUS is not homogeneous). However, at this point in the whole
# selection process, we don't have enough sources in order to properly compute magnitude cuts corresponding to a certain
# S/N threshold for each tile. Instead of looking for the S/N-dependant mag_limit cut, then, I am directly comparing
# the S/N of each source to a certain threshold (S/N=3) in order to include/exclude them from the final list.

SNcutR_mask = ( (1./jpl['rJAVA'][:,1]) > setup['SN'])
jpl['over_rJAVA_cut'] = SNcutR_mask

SNcutG_mask = ( (1./jpl['gJAVA'][:,1]) > setup['SN'])
jpl['over_gJAVA_cut'] = SNcutG_mask

jpl['over_all_SNcuts'] = ((SNcutR_mask) & (SNcutG_mask))




'''
#-------------- RAFA CATALOG CROSSMATCH  --------------#
if setup['load_rafa'] == True:
    dist_rafa, ind_rafa = xmatch(jpl['coords'], rafa['coords'], max_distance=3./3600.)
    rafa_mask = (dist_rafa != np.inf)
    print 'HII regions: ', len(jpl['coords'][rafa_mask,0])
else:
    rafa_mask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
jpl['HII_region'] = rafa_mask

hII_extended_mask = ((jpl['extended']) & (jpl['HII_region']))
print 'There are ', len(jpl['coords'][hII_extended_mask,0]),' HII regions in the list of extended candidates'
'''



#-------------- CDS CROSSMATCH SERVICE (GALEX) --------------#
if setup['galexmask'] == True:
    dist_gj, ind_gj = xmatch(jpl['coords'], galex['coords'], max_distance=3./3600.)
    galex_mask= (dist_gj != np.inf)
    # cc = 0
    # for i in np.arange(len(jpl['coords'][:,0])):
    #     tools.plot_JPLUS_results.plot_JPLUSphotoSpectra(jpl, i, mask=[], units='mags',\
    #                                                     zfromSDSS=False, number=cc)
    #     cc += 1
    print 'Not in galex x-match: ', len(jpl['coords'][~galex_mask,0])
else:
    galex_mask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
jpl['in_galex'] = galex_mask




#-------------- SDSS GALAXIES CROSSMATCH  --------------#
if setup['load_sdssGal'] == True:
    dist_gal, ind_gal = xmatch(jpl['coords'], sdss_gal['coords'], max_distance=3./3600.)
    sdss_gal_mask = (dist_gal != np.inf)
    sdss_gal_num =  len(jpl['coords'][sdss_gal_mask,0])
    sdss_gal_contam = int(10000.*sdss_gal_num/( len(jpl['coords'][:,0])) )/100.
    
    jpl['redshift'][sdss_gal_mask] = sdss_gal['zspec'][ind_gal[sdss_gal_mask]]
    print 'Number of sdss galaxies (all z): ', sdss_gal_num
    print 'sdss_galaxies contamination: %2.2f'%sdss_gal_contam+'% (all-z galaxies)'

    #---- ANALYSIS on SDSS GALAXIES ONLY ----#
    # sdss_gal_inJPL = jplus.tools.select_object(jpl, sdss_gal_mask)
    # filts = ['J0378', 'uJAVA', 'gJAVA']
    # J0378_mask_sdss_gal = tools.find_lineExcess(sdss_gal_inJPL, fset=filts)
    # filts = ['J0660', 'rJAVA', 'iJAVA']
    # J0660_mask_sdss_gal = tools.find_lineExcess(sdss_gal_inJPL, fset=filts)
    # sdss_gal_linemask = ((J0378_mask_sdss_gal) & (J0660_mask_sdss_gal))
    # print 'Number of sdss galaxies emitting in J0378+J0660: ',\
    #     len(sdss_gal_inJPL['coords'][sdss_gal_linemask,0])
    # # plt.hist(sdss_gal_inJPL['redshift'], color='b', alpha=0.5, label='all sdss gal in JPLUS')
    # # plt.hist(sdss_gal_inJPL['redshift'][sdss_gal_linemask], color='r', label='sdss liners')
    # # plt.show()
    # # plt.close()
    # #from tools import plot_JPLUSphotoSpectra as JPLphotSP
    # #JPLphotSP(sdss_gal_inJPL, 0, mask=[], units='flux', zsdss=0.36, zfromSDSS=True, number=0)
else:
    sdss_gal_mask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
jpl['sdss_gal'] = sdss_gal_mask




#-------------- SDSS STARS CROSSMATCH  --------------#
if setup['load_sdssSTAR'] == True:
    dist_str, ind_str = xmatch(jpl['coords'], sdss_str['coords'], max_distance=3./3600.)   # xmatch with candidates
    starmask = (dist_str != np.inf)
    starnum =  len(jpl['coords'][starmask,0])
    nstr = 100.*starnum/(len(jpl['coords'][:,0]))
    
    print 'Number of sdss stars: ', starnum
    print 'sdss_stars contamination: %2.2f'%nstr+'% '
else:
    starmask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
jpl['sdss_stars'] = starmask




#-------------- SDSS QUASARS CROSSMATCH  --------------#
if setup['load_sdssQSO'] == True:
    # dist, ind = xmatch(alljpl['coords'],qso['coords'],max_distance=3./3600.) # xmatch with all jplus
    dist, ind = xmatch(jpl['coords'], qso['coords'], max_distance=3./3600.)   # xmatch with candidates
    qsomask = (dist != np.inf)
    qsonum =  len(jpl['coords'][qsomask,0])
    n = 100.*qsonum/(len(jpl['coords'][:,0]))
    jpl['redshift'][qsomask] = qso['zspec'][ind[qsomask]]   #using sdss redshift for jplus crossmatched sources

    dist_z, ind_z = xmatch(jpl['coords'], qso_z['coords'], max_distance=3./3600.) # right-z quasars
    qsomask_z = (dist_z != np.inf)
    qsonum_z =  len(jpl['coords'][qsomask_z,0])
    n_z = 100.*qsonum_z/(len(jpl['coords'][:,0]))

    print 'Number of sdss quasars: ', qsonum
    print 'Number of sdss quasars at the right z: ', qsonum_z
    print 'sdss_quasars contamination: %2.2f'%n+'% (all-z quasars)'

    # tools.plot_quasar_ColorColor(jpl, qso, qsomask, ind)
    # tools.plot_quasar_zDistrib(qso, qsomask, ind)
    # tools.plot_quasar_PhotoSpec(jpl, qso, qsomask, ind)
    # tools.plot_quasar_lumifunc(jpl, qso, qsomask, ind, n_tiles=149)
    # # BE CAREFUL ABOUT n_tiles (it should be the TOTAL number of OBSERVED tiles in the jplus survey)
else:
    qsomask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
    qsomask_z = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
    
jpl['sdss_qso'] = qsomask
jpl['sdss_rightZ_qso'] = qsomask_z




#-------------- LQAC CATALOG CROSSMATCH  --------------#
if setup['load_lqac'] == True:
    dist_lqac, ind_lqac = xmatch(jpl['coords'], lqac['coords'], max_distance=3./3600.)
    lqac_mask = (dist_lqac != np.inf)
else:
    lqac_mask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
jpl['lqac_qso'] = lqac_mask

only_lqac_mask = ((jpl['lqac_qso']) & (~jpl['sdss_qso']))
print 'New quasars from LQAC (not already in sdss): ', len(jpl['coords'][only_lqac_mask,0])




#-------------- GAIA CROSSMATCH  --------------#
if setup['load_gaia']==True:
    djg, ijg = xmatch(jpl['coords'], gaia['coords'], max_distance=3./3600.)
    gaiamask = (djg != np.inf)
    
    gaianum = len(jpl['coords'][gaiamask,0])
    ngaia = 100.*gaianum/(len(jpl['coords'][:,0]))

    print 'Number of gaia objects: ', gaianum
    print 'gaia contamination: %2.2f'%ngaia+'% (all gaia objs)'
else:
    gaiamask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
jpl['gaia'] = gaiamask




# #-------------- CANDIDATE EXTENT  --------------#
print 'Extended objects: ', len(jpl['coords'][jpl['extended'],0])
print 'Compact objects: ', len(jpl['coords'][jpl['compact'],0])




#--------------  (QSO?) NB LINE EMISSIONS   --------------#
#tools.plot_JPLUSfilters_QSOlines()

#-- J0660 --#
filts = ['J0660', 'rJAVA', 'iJAVA']
J0660_mask = tools.find_lineExcess(jpl, fset=filts)
jpl['J0660_liners'] = J0660_mask

#-- J0861 --#
filts = ['J0861', 'iJAVA', 'zJAVA']
J0861_mask = tools.find_lineExcess(jpl, fset=filts)
jpl['J0861_liners'] = J0861_mask

#-- J0410 --#
filts = ['J0410', 'uJAVA', 'gJAVA']
J0410_mask = tools.find_lineExcess(jpl, fset=filts)
jpl['J0410_liners'] = J0410_mask

#-- J0515 --#
filts = ['J0515', 'gJAVA', 'rJAVA']
J0515_mask = tools.find_lineExcess(jpl, fset=filts)
jpl['J0515_liners'] = J0515_mask

#-- J0430 --#
filts = ['J0430', 'uJAVA', 'gJAVA']
J0430_mask = tools.find_lineExcess(jpl, fset=filts)
jpl['J0430_liners'] = J0430_mask

lines_mask = ( (J0410_mask) | (J0430_mask) | (J0515_mask) | (J0660_mask) | (J0861_mask) )
print 'Line emitters: ', len(jpl['coords'][lines_mask,0])

jpl['all_liners'] = lines_mask



#----------------------  Lya LUMINOSITIES  ----------------------#
fluxlya, fcont, flux_nb, fcont_error = tools.threeFM.jplus3FM(jpl)
comov_table = setup['tools'] + 'comoving_distance.txt'
z, eta = np.genfromtxt(comov_table, skip_header=1, unpack=True)
line = interp1d(z, eta, kind='nearest', bounds_error=None, fill_value='extrapolate')
dl = line(setup['zmean'])*(1+setup['zmean'])                   #now in Mpc, all at same redshift
dl = dl*(1.e6)*3.0856e18              #now in cm 
lumilya = 4.*np.pi*fluxlya*(dl**2.)   #now in erg/s

jpl['lya_lumin'] = lumilya
print min(lumilya), max(lumilya)




#------------------  CANDIDATES SELECTION  ------------------#
candidates = ( (~jpl['in_galex']) & (~jpl['sdss_gal']) & (~jpl['sdss_qso']) & (~jpl['sdss_stars']) & \
              (~jpl['lqac_qso']) & (~jpl['all_liners']) & (~jpl['gaia']) )

if setup['morph_sel'] == 'extd':
    candidates = ( (candidates) & (jpl['extended']) )
else:
    candidates = ( (candidates) & (jpl['compact']) )

    
if setup['morph_sel'] == 'extd': print '\n---------- EXTENDED ----------'
if setup['morph_sel'] == 'comp': print '\n---------- COMPACT ----------'
print 'Candidates are defined as:'
if setup['galexmask'] == True: print '- NOT in Galex'
if setup['load_sdssGal'] == True: print '- NOT sdss galaxies'
if setup['load_sdssQSO'] == True: print '- NOT sdss QSOs'
if setup['load_sdssSTAR'] == True: print '- NOT sdss stars'
if setup['load_lqac'] == True: print '- NOT in LQAC'
if setup['load_gaia'] == True: print '- NOT gaia stars'
print '- NOT "liners"'
if setup['morph_sel'] == 'extd': print 'EXTENDED objects\n'
if setup['morph_sel'] == 'comp': print 'COMPACT objects\n'

over_SN = ((candidates) & (jpl['over_all_SNcuts']))
print 'Final candidates number: ', len(jpl['coords'][candidates,0]),'; (', len(jpl['coords'][over_SN,0]),'of which are above all BB S/N cuts)','\n'




#--------------  LOW-Z OII EMITTERS   --------------#
if setup['filters'][0] == 'J0378':
    
    #------- COLOR-COLOR PLOT -------#
    clrx = ['gJAVA', 'rJAVA']
    clry = ['iJAVA', 'zJAVA']
    OIIx = jpl[clrx[0]][jpl['J0660_liners'],0] - jpl[clrx[1]][jpl['J0660_liners'],0]
    OIIy = jpl[clry[0]][jpl['J0660_liners'],0] - jpl[clry[1]][jpl['J0660_liners'],0]
    allx = jpl[clrx[0]][:,0] - jpl[clrx[1]][:,0]
    ally = jpl[clry[0]][:,0] - jpl[clry[1]][:,0]
    candsx = jpl[clrx[0]][candidates,0] - jpl[clrx[1]][candidates,0]
    candsy = jpl[clry[0]][candidates,0] - jpl[clry[1]][candidates,0]
    f = plt.figure(figsize=(12,10))
    plt.plot(OIIx, OIIy, 'ob', markersize=6, alpha=0.8, label='OII_em')
    plt.plot(candsx, candsy, 'or', markersize=6, alpha=0.4, label='candidates')
    #plt.plot(allx, ally, 'og', markersize=3, alpha=0.2, label='all')
    plt.xlabel(clrx[0]+' - '+clrx[1]+'  [mags]', fontsize=14)
    plt.ylabel(clry[0]+' - '+clry[1]+'  [mags]', fontsize=14)
    plt.legend(fontsize=14, fancybox=True)
    plt.title('OII vs candidates  color plot', fontsize=17)
    #plt.savefig(setup['plots']+'OII_emitters/OIIemitrs_'+clrx[0][0]+clrx[1][0]+'-'+clry[0][0]+clry[1][0]+'JAVA.png')
    #plt.show()
    plt.close()

    #------- COMPOSITE SPECTRUM -------#
    tools.plot_JPLUS_results.plot_MEDspectra(jpl, mask=jpl['J0660_liners'], titol='OII emitters median photo-spectrum')
    tools.plot_JPLUS_results.plot_MEDspectra(jpl, mask=jpl['extended'], titol='extended candidates median photo-spectrum')
    
    #------- SOME STATISTICS -------#
    print '------- OII emittters -------'
    print 'Total number: ',len(jpl['rJAVA'][jpl['J0660_liners'],0])

    test = ((jpl['J0660_liners']) & (jpl['extended']))
    print 'Extended: ',len(jpl['rJAVA'][test,0])
    test = ((jpl['J0660_liners']) & (jpl['compact']))
    print 'Compact: ',len(jpl['rJAVA'][test,0])

    #frac = float( len(jpl['rJAVA'][jpl['J0660_liners'],0]) )/len(jpl['rJAVA'][jpl['extended'],0])
    #print 'Fraction over extended: %2.2f' %frac
    frac = float(len(jpl['rJAVA'][jpl['J0660_liners'],0]))/len(jpl['rJAVA'][:,0])
    print 'Fraction over all: %2.2f'%frac

    #------- GALEX CROSSMATCH -------#
    matrix = np.empty( (2, len(jpl['coords'][jpl['J0660_liners'],0])) )      # RA-DEC list of OII emitters
    matrix[0,:] = jpl['coords'][jpl['J0660_liners'],0]
    matrix[1,:] = jpl['coords'][jpl['J0660_liners'],1]
    fmt_string = '%3.5f %3.5f'
    #np.savetxt(setup['results']+'OII-emitters_RA-DEC_list.csv', matrix.T, fmt=fmt_string)
    gal_ra,gal_dec = np.genfromtxt(setup['data']+'OII-emitters_GALEX-(All_IR_Sky)_xmatch_19-06-2017.csv', delimiter=',', skip_header=1, usecols=(1,2), unpack=True)
    print 'In galex: '+str(1)+' over '+str(len(jpl['rJAVA'][jpl['J0660_liners'],0]))+'\nxmatch coords: ',str(gal_ra.T)+' '+str(gal_dec.T)
else:
    '\nNo OII_emitters analysis, we need J0378 filter!\n'




    
#--------   GET THE QUASARS BY THEIR COLORS!   --------#
print '\n\n'
#add_mask = ((jpl['extended']) & (~jpl['sdss_qso']) )
com_mask = ((jpl['compact']) & (~jpl['sdss_qso']) & (~jpl['all_liners']) & (~jpl['in_galex']) & (~jpl['sdss_gal']) & (~jpl['sdss_stars']) & (~jpl['lqac_qso']) & (~jpl['all_liners']) & (~jpl['gaia']) )
ext_mask = ((jpl['extended']) & (~jpl['sdss_qso']) & (~jpl['all_liners']) & (~jpl['in_galex']) & (~jpl['sdss_gal']) & (~jpl['sdss_stars']) & (~jpl['lqac_qso']) & (~jpl['all_liners']) & (~jpl['gaia']) )
plt.figure(figsize=(12,10))
fx = ['gJAVA', 'rJAVA']
fy = ['gJAVA', 'iJAVA']
xcol = jpl[fx[0]][:,0] - jpl[fx[1]][:,0]
ycol = jpl[fy[0]][:,0] - jpl[fy[1]][:,0]
#--- color-cut criterion ---#
on_QSO = True  #if 'True' the following color cut will be fine-tuned in order to select SDSS QSOs with high completeness and purity, if 'False' the color cut will be defined on extended emitters
if (setup['filters'][0] == 'J0378') or (setup['filters'][0] == 'J0395'):
    if on_QSO == True:
        for i in np.arange(0, 1.6, 0.2):
            coul = (i/1.41, 0.1, 0.8) #, 0.7, 0.4)
            cut = ((xcol < i) & (ycol < i))    # QSOs criterion
            qsos_above = ((cut) & (jpl['sdss_qso']) )
            sources_above = ((cut) & (ext_mask))
            complet = round( float(len(xcol[qsos_above]))/len(xcol[jpl['sdss_qso']])*100. , 2)
            purit = round( float(len(xcol[qsos_above]))/(len(xcol[qsos_above]) + len(xcol[sources_above]))*100. , 2)
            plt.plot((i, i), (-0.9, i), color=coul, linestyle='--', linewidth=1.5)
            plt.plot((-2.4, i), (i, i), color=coul, linestyle='--', linewidth=1.5)
            plt.text(i-0.08, -1., 'cut:', color=coul, fontweight='bold', fontsize=8)
            plt.text(i-0.08, -1.1, str(i), color=coul, fontweight='bold', fontsize=8)
            plt.text(-2.4, i+0.03, 'Comp: '+str(complet)+'%', color=coul, fontweight='bold', fontsize=9)
            plt.text(-1.5, i+0.03, 'Pur: '+str(purit)+'%', color=coul, fontweight='bold', fontsize=9)
        # another limit, by hand
        cut = ((xcol < 0.5) & (ycol < 1.0))    # QSOs criterion
        qsos_above = ((cut) & (jpl['sdss_qso']) )
        sources_above = ((cut) & (ext_mask))
        complet = round( float(len(xcol[qsos_above]))/len(xcol[jpl['sdss_qso']])*100. , 2)
        purit = round( float(len(xcol[qsos_above]))/(len(xcol[qsos_above]) + len(xcol[sources_above]))*100. , 2)
        plt.plot((0.5, 0.5), (-1.5, 1.0), 'k--', linewidth=1.5)
        plt.plot((-2.4, 0.5), (1.0, 1.0), 'k--', linewidth=1.5)
        plt.text(0.55, -1.5, 'Comp: '+str(complet)+'%'+'  Pur: '+str(purit)+'%', color='k', fontweight='bold', fontsize=10)
    else:
        for i in np.arange(0, 1.1, 0.2):
            coul = (i, 0.7, 0.4)
            cut = ((xcol > i) & (ycol > i))
            qsos_above = ((cut) & (jpl['sdss_qso']) )
            sources_above = ((cut) & (ext_mask))
            complet = round( float(len(xcol[sources_above]))/len(xcol[ext_mask])*100. , 2)
            purit = round( float(len(xcol[sources_above]))/(len(xcol[qsos_above]) + len(xcol[sources_above]))*100. , 2)
            plt.plot((i, i), (i, 3.5), color=coul, linestyle='--', linewidth=2)
            plt.plot((i, 2.6), (i, i), color=coul, linestyle='--', linewidth=2)
            plt.text(i, 3.7, 'cut:', color=coul, fontweight='bold', fontsize=9)
            plt.text(i, 3.6, str(i), color=coul, fontweight='bold', fontsize=9)
            plt.text(2.7, i, 'Comp: '+str(complet)+  '   Pur: '+str(purit), color=coul, fontweight='bold', fontsize=9)
        plt.plot(xcol[~jpl['sdss_qso']], ycol[~jpl['sdss_qso']], 'og', markersize=2, alpha=0.3, label='all NOT-QSOs')  # ALL not-QSO sources
        
plt.plot(xcol[jpl['sdss_rightZ_qso']], ycol[jpl['sdss_rightZ_qso']], 'or', markersize=4, alpha=0.6, label='SDSS right_z QSOs')
not_rightZ_qso = ((jpl['sdss_qso']) & (~jpl['sdss_rightZ_qso']))
plt.plot(xcol[not_rightZ_qso], ycol[not_rightZ_qso], 'xk', markersize=4, alpha=0.98, label='SDSS "wrong"_z QSOs')
plt.plot(xcol[ext_mask], ycol[ext_mask], 'ob', markersize=4, alpha=0.8, label='extended candidates')
plt.plot(xcol[com_mask], ycol[com_mask], 'og', markersize=2, alpha=0.6, label='compact candidates')
plt.xlabel(fx[0]+' - '+fx[1]+'  [mags]')
plt.ylabel(fy[0]+' - '+fy[1]+'  [mags]')
plt.legend(fancybox=True)
plt.title('Color-color cut to discard QSOs (data from '+setup['filters'][0]+' filter)', fontsize=16)
# if on_QSO == True:
#     plt.savefig(setup['plots']+fx[0][0]+fx[1][0]+'-'+fy[0][0]+fy[1][0]+'JAVA_QSOs-vs-extended_ColorColor_withCutsOnQSOs_withCompact'+setup['filters'][0]+'.png')
#     a = 1
# else:
#     plt.savefig(setup['plots']+fx[0][0]+fx[1][0]+'-'+fy[0][0]+fy[1][0]+'JAVA_QSOs-vs-extended_ColorColor_withCuts_'+setup['filters'][0]+'.png')
#     a = 1
#plt.show()
plt.close()
#sys.exit()







#-------------------------   PLOTS   -------------------------#
# #FLAG-ORIENTED COLOR-COLOR
# fig = plt.figure(figsize=(12,10))
# fx = ['gJAVA', 'rJAVA']
# fy = ['rJAVA', 'iJAVA']
# xcol = jpl[fx[0]][:,0] - jpl[fx[1]][:,0]
# ycol = jpl[fy[0]][:,0] - jpl[fy[1]][:,0]
# plt.plot(xcol[jpl['sdss_qso']], ycol[jpl['sdss_qso']], 'or', markersize=4, alpha=0.3, label='all SDSS QSOs')
# #plt.plot(xcol[jpl['sdss_rightZ_qso']], ycol[jpl['sdss_rightZ_qso']], 'om', markersize=4, alpha=0.8, label='SDSS QSOs z=2.1')
# plt.plot(xcol[~jpl['sdss_qso']], ycol[~jpl['sdss_qso']], 'og', markersize=2, alpha=0.3, label='all NOT-QSOs')
# plt.xlabel(fx[0]+' - '+fx[1]+'  [mags]')
# plt.ylabel(fy[0]+' - '+fy[1]+'  [mags]')
# #plt.axis([-1, 2, -0.5, 2.5])
# plt.legend(fancybox=True)
# plt.savefig(setup['plots']+fx[0][0]+fx[1][0]+'-'+fy[0][0]+fy[1][0]+'JAVA_QSOs-vs-extended_ColorColor.png')
# plt.show()
# plt.close()
# sys.exit()

# #PLOT SINGLE SPECTRA
# objectid = 6
# tools.plot_JPLUS_results.plot_JPLUSphotoSpectra(jpl, objectid-1, mask=jpl['sdss_rightZ_qso'], units='flux', zfromSDSS=False)

# #PLOT MEDIAN SPECTRA ("""COMPOSITE""")
# mask_med_spec = jpl['sdss_rightZ_qso']
# tools.plot_JPLUS_results.plot_MEDspectra(jpl, mask=mask_med_spec, titol='QSOs SDSS median Photo-spectrum')

# candidates2 = ( (~jpl['in_galex']) & (~jpl['sdss_gal']) & (~jpl['sdss_qso']) & \
#                (~jpl['sdss_rightZ_qso']) & (jpl['extended']) & (~jpl['all_liners']) &\
#                (~jpl['gaia']) )
# mask_med_spec2 = candidates2

# tools.plot_JPLUS_results.plot_MEDspectra(jpl, mask=mask_med_spec2, titol='Extended candidates median Photo-spectrum')
# tools.plot_JPLUS_results.plot_MEDspectra_diff(jpl, mask1=mask_med_spec, mask2=mask_med_spec2, titol='Difference median Photo-spectrum')





#--------------------------------  DATA SAVING  --------------------------------#
all_final = jplus.tools.select_object(jpl, candidates)       # all candidates (excluding QSOs, galaxies, liners... see definition of "candidates" mask, above)
final = jplus.tools.select_object(jpl, over_SN)              # same as above but including only candidates above Broad-Band S/N cuts

#HDF5 catalogues
dd.io.save(setup['final_catalog'], all_final)                # all candidates (excluding QSOs, galaxies, liners...)
dd.io.save(setup['final_catalog'][:-3]+'_overSN.h5', final)  # only candidates above Broad-Band S/N cuts
dd.io.save(setup['flagged_catalog'], jpl)                    # ALL candidates (excluding nothing: this rewrites the input, adding the selection flags)

#ASCII tables
matrix = np.empty( (8, len(jpl['coords'][candidates,0])) )   # all candidates (excluding QSOs, galaxies, liners...)
matrix[0,:] = jpl['tile_id'][candidates]
matrix[1,:] = jpl['coords'][candidates,0]
matrix[2,:] = jpl['coords'][candidates,1]
matrix[3,:] = jpl[setup['filters'][0]][candidates,0]
matrix[4,:] = fluxlya[candidates]
matrix[5,:] = lumilya[candidates]
matrix[6,:] = jpl['rJAVA'][candidates,0]
matrix[7,:] = jpl['redshift'][candidates]
head = 'tl   ra        dec    NBmag line_flux line_lum  rJAVA  z_NB'
fmt_string = '%4d %3.5f %3.5f %4.2f %4.3e %4.3e %4.2f %4.3f'
np.savetxt(setup['final_list'], matrix.T, fmt=fmt_string, header=head)

matrix = np.empty( (2, len(jpl['coords'][candidates,0])) )   # RA-DEC of all candidates (excluding QSOs, galaxies, liners...)
matrix[0,:] = jpl['coords'][candidates,0]
matrix[1,:] = jpl['coords'][candidates,1]
head = 'ra        dec'
fmt_string = '%3.5f %3.5f'
np.savetxt(setup['final_radec'], matrix.T, fmt=fmt_string, header=head)

matrix = np.empty( (2, len(jpl['coords'][over_SN,0])) )      # RA-DEC of candidates above Broad-Band S/N cuts
matrix[0,:] = jpl['coords'][over_SN,0]
matrix[1,:] = jpl['coords'][over_SN,1]
head = 'ra        dec'
fmt_string = '%3.5f %3.5f'
np.savetxt(setup['final_radec_overSN'], matrix.T, fmt=fmt_string, header=head)

sys.exit()


#--------------- PLOT LUMINOSITY FUNCTION ---------------#
n_tiles = 205
jplus_area = (1.4**2.)*n_tiles #now in sq.degrees
sphere_area = 4*np.pi/((np.pi/180.)**2.) #now in sq.degrees
fsky = jplus_area/sphere_area

propVol = False
tools.plot_JPLUS_results.plot_lumiFunc(fsky, line, propVol)

#sys.exit()



# FOR ALESSANDRO:
# jpl_foot_mask = ((sdss_gal['coords'][:,0] > 126.25) & (sdss_gal['coords'][:,0] < 141.25) &\
#                  (sdss_gal['coords'][:,1] > 51.) & (sdss_gal['coords'][:,1] < 60.))
# sdss_for_ale = jplus.tools.select_object(sdss_gal, jpl_foot_mask)

# dist_sg, ind_sg = xmatch(sdss_for_ale['coords'], gaia['coords'], max_distance=3./3600.)
# mask_sg = (dist_sg != np.inf)
# print 'xmatch: ', len(sdss_for_ale['coords'][mask_sg,0])
# print 'gaia: ', len(gaia['coords'][:,0])
# print 'sdss: ', len(sdss_for_ale['coords'][:,0])
