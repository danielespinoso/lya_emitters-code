import sys
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
jpl = dd.io.load(setup['jplus_candidates'])  # loads the list of candidates
#alljpl = dd.io.load(setup['jplus_input'])   # loads the initial jplus catalogue
    
if setup['load_sdssGal'] == True:
    sdss_gal = dd.io.load(setup['sdss_gal'])  #load sdss galaxies catalogue

if setup['load_sdssQSO'] == True:
    qso = dd.io.load(setup['sdss_qso'])
    qso_zmask = ((qso['redshift'] > setup['zmin']) & (qso['redshift'] < setup['zmax']))
    qso_z = jplus.tools.select_object(qso, qso_zmask)
    
if setup['load_gaia'] == True:
    gaia = dd.io.load(setup['gaia'])          #load gaia stars catalogue

if setup['load_mock_cd'] == True:
    mock = dd.io.load(setup['mock_candidates'])
    
if setup['galexmask'] == True:
    galex = dd.io.load(setup['galex_in'])
    
if setup['load_rafa'] == True:
    rafa = dd.io.load(setup['rafa'])          #load rafa catalog of HII regions in jplus
    
if setup['load_lqac'] == True:
    lqac = dd.io.load(setup['lqac'])          #load Large Qasar Astrometric Catalogue

print 'Loaded datasets:\n'
if setup['load_sdssGal'] == True: print 'sdss galaxies'
if setup['load_sdssQSO'] == True: print 'sdss quasars'
if setup['load_gaia'] == True: print 'overlapping gaia'
if setup['load_mock_cd'] == True: print 'jplus mocks'
if setup['galexmask'] == True: print 'galex sources'
if setup['load_rafa'] == True: print 'rafa HII sources'
if setup['load_lqac'] == True: print 'lqac catalogue'

print '\nInitial candidates: ', len(jpl['coords'][:,0])


# FOR ALESSANDRO:
# jpl_foot_mask = ((sdss_gal['coords'][:,0] > 126.25) & (sdss_gal['coords'][:,0] < 141.25) &\
#                  (sdss_gal['coords'][:,1] > 51.) & (sdss_gal['coords'][:,1] < 60.))
# sdss_for_ale = jplus.tools.select_object(sdss_gal, jpl_foot_mask)

# dist_sg, ind_sg = xmatch(sdss_for_ale['coords'], gaia['coords'], max_distance=3./3600.)
# mask_sg = (dist_sg != np.inf)
# print 'xmatch: ', len(sdss_for_ale['coords'][mask_sg,0])
# print 'gaia: ', len(gaia['coords'][:,0])
# print 'sdss: ', len(sdss_for_ale['coords'][:,0])




#-------------- BB SIGNAL-to-NOISE SELECTION  --------------#
mcutR = tools.jplus_SNratio(jpl, mask=[], band='rJAVA', SNcut=3.)
mcutU = tools.jplus_SNratio(jpl, mask=[], band='uJAVA', SNcut=3.) #something strange happens
mcutU = 20.15
mcutG = tools.jplus_SNratio(jpl, mask=[], band='gJAVA', SNcut=3.)

SNcutR_mask = (jpl['rJAVA'][:,0] < mcutR)
jpl['over_rJAVA_cut'] = SNcutR_mask

SNcutU_mask = (jpl['uJAVA'][:,0] < mcutU)
jpl['over_uJAVA_cut'] = SNcutU_mask

SNcutG_mask = (jpl['gJAVA'][:,0] < mcutG)
jpl['over_gJAVA_cut'] = SNcutG_mask

jpl['over_all_SNcuts'] = ((SNcutR_mask) & (SNcutU_mask) & (SNcutG_mask))


#---------------  REDSHIFT ESTIMATE and SELECTION  --------------------#
#------- LEPHARE -------#
#jpl_LPz = lya_lephare_jplus(jpl, setup['leph_bb'], name=' ', owtspec=True, owr_prep=True, owr_run=True)
# leph_results = setup['data'] + 'jplus_alltoR24_autoMags_upad_dual_jplusPhot_lpphotoz_recalpertile_filflag0000010101010000000000_EM_ .h5'
# leph = dd.io.load(leph_results)

# plt.plot(leph['PDFzstep'], leph['PDF'][4,:], 'r-')
# plt.xlabel('z')
# plt.ylabel('PDF')
# plt.show()
# plt.close()

'''
jpl['redshift'] = zml_NB  ### USE LEPHARE ESTIMATE
jpl['redshift_noNB'] = zml  ### USE LEPHARE ESTIMATE
'''
if setup['zmask'] == True:
    zmask = ((jpl['redshift'] > setup['zmin']) & (jpl['redshift'] < setup['zmax']))
else:
    zmask = np.ones(len(jpl['coords'][:,0]), dtype=bool)
jpl['right_z'] = zmask

print 'After redshift selection: ', len(jpl['coords'][zmask,0])




#-------------- RAFA CATALOG CROSSMATCH  --------------# (TRY)
if setup['load_rafa'] == True:
    dist_rafa, ind_rafa = xmatch(jpl['coords'], rafa['coords'], max_distance=3./3600.)
    rafa_mask = (dist_rafa != np.inf)
    print 'HII regions: ', len(jpl['coords'][rafa_mask,0])
else:
    rafa_mask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
jpl['HII_region'] = rafa_mask




#-------------- LQAC CATALOG CROSSMATCH  --------------# (TRY)
if setup['load_lqac'] == True:
    dist_lqac, ind_lqac = xmatch(jpl['coords'], lqac['coords'], max_distance=3./3600.)
    lqac_mask = (dist_lqac != np.inf)
    print 'LQAC quasars: ', len(jpl['coords'][lqac_mask,0])
else:
    lqac_mask = np.zeros(len(jpl['coords'][:,0]), dtype=bool)
jpl['lqac_qso'] = lqac_mask




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
    
    jpl['redshift'][sdss_gal_mask] = sdss_gal['redshift'][ind_gal[sdss_gal_mask]]
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




#-------------- SDSS QUASARS CROSSMATCH  --------------#
if setup['load_sdssQSO'] == True:
    # dist, ind = xmatch(alljpl['coords'],qso['coords'],max_distance=3./3600.) # xmatch with all jplus
    dist, ind = xmatch(jpl['coords'], qso['coords'], max_distance=3./3600.)   # xmatch with candidates
    qsomask = (dist != np.inf)
    qsonum =  len(jpl['coords'][qsomask,0])
    n = 100.*qsonum/(len(jpl['coords'][:,0]))

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




#-------------- CANDIDATE EXTENT  --------------# TO BE UPDATED AFTER THE NEW MORPHOLOGY SELECTION
pxl = 0.55*0.55  # jplus pixel area in arcsec^2
mumax = jpl['mu_max_r'] - 2.5*np.log10(pxl)
morph = mumax - jpl['rJAVA'][:,0]

tile_ext_lim = np.genfromtxt('tile_to_extent.txt', unpack=False)
jpl['compactness'] = np.zeros(len(jpl['mu_max_r']), dtype=bool)
for i in set(jpl['tile_id']):
    ext_lim = tile_ext_lim[tile_ext_lim[:,0]==i, 1][0]
    mask_tile = (jpl['tile_id'] == i)
    mask_tilext = ((mask_tile) & (morph < ext_lim))
    
    jpl['compactness'] = ((jpl['compactness']) | (mask_tilext))  #EXTENT FLAG (True if compact)

print 'Extended objects: ', len(jpl['coords'][~jpl['compactness'],0])




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






#------------------  FLAGS AND CANDIDATE SELECTION ------------------#
flags = ['right_z', 'in_galex', 'sdss_gal', 'sdss_qso', 'sdss_rightZ_qso', 'compactness',\
         'J0660_liners', 'J0861_liners', 'J0410_liners', 'J0430_liners', 'J0515_liners',\
         'all_liners', 'gaia', 'HII_region', 'lqac_qso']
candidates = ( (jpl['right_z']) & (~jpl['in_galex']) & (~jpl['sdss_gal']) & (~jpl['sdss_qso']) & \
               (~jpl['sdss_rightZ_qso']) & (~jpl['compactness']) & (~jpl['all_liners']) &\
               (~jpl['gaia']) )

print 'Final candidates number: ', len(jpl['coords'][candidates,0]),'\n\n'

over_SN = ((candidates) & (jpl['over_all_SNcuts']))
print 'Over all SN cuts:', len(jpl['coords'][over_SN,0]),'\n\n'

sdss_lqac_mask = ((jpl['lqac_qso']) & (jpl['sdss_qso']))
print 'There are ', len(jpl['coords'][sdss_lqac_mask,0]),' lqac quasars in sdss list'
hII_extended_mask = ((~jpl['compactness']) & (jpl['HII_region']))
print 'There are ', len(jpl['coords'][hII_extended_mask,0]),' HII regions in the extended list'

# fig = plt.figure(figsize=(12,10))
# plt.plot(jpl['rJAVA'][jpl['sdss_gal'],0], morph[jpl['sdss_gal']], 'oc', markersize=4, label='SDSS gal', alpha=1.0)
# plt.plot(jpl['rJAVA'][jpl['sdss_qso'],0], morph[jpl['sdss_qso']], 'ob', markersize=4, label='SDSS QSO (all)', alpha=0.5)
# plt.plot(jpl['rJAVA'][jpl['sdss_rightZ_qso'],0], morph[jpl['sdss_rightZ_qso']], 'or', markersize=4,label='SDSS QSO (z=2.1)', alpha=1.)
# plt.plot(jpl['rJAVA'][candidates,0], morph[candidates], 'og', markersize=2, label='candidates', alpha=0.8)
# plt.legend()
# plt.xlabel('rJAVA  [mags]')
# plt.ylabel('MU_MAX(rJAVA)  -  rJAVA'+'  [mag]')
# plt.title('All JPLUS tiles')
# plt.show()
# plt.close()


# #------- FLAG-ORIENTED COLOR-COLOR PLOT -------#
# fig = plt.figure(figsize=(12,10))
# fx = ['uJAVA', 'gJAVA']
# fy = ['gJAVA', 'rJAVA']
# xcol = jpl[fx[0]][:,0] - jpl[fx[1]][:,0]
# ycol = jpl[fy[0]][:,0] - jpl[fy[1]][:,0]
# #plt.plot(xcol[jpl['sdss_qso']], ycol[jpl['sdss_qso']], 'or', markersize=4, alpha=0.5, label='all SDSS QSOs')
# plt.plot(xcol[jpl['sdss_rightZ_qso']], ycol[jpl['sdss_rightZ_qso']], 'om', markersize=4, alpha=0.8, label='SDSS QSOs z=2.1')
# plt.plot(xcol[over_SN], ycol[over_SN], 'ob', markersize=4, alpha=0.7, label='extd over all SNcuts')
# #plt.plot(xcol[jpl['gaia']], ycol[jpl['gaia']], 'og', markersize=4, alpha=0.5, label='gaia')
# plt.xlabel(fx[0]+' - '+fx[1]+'  [mags]')
# plt.ylabel(fy[0]+' - '+fy[1]+'  [mags]')
# #plt.axis([-1, 2, -0.5, 2.5])
# plt.legend()
# plt.show()
# plt.close()

#sys.exit()




#--------------  PLOT SPECTRA (composite is WIP)  --------------#
# objectid = 6
# tools.plot_JPLUS_results.plot_JPLUSphotoSpectra(jpl, objectid-1, mask=jpl['sdss_rightZ_qso'], units='flux', zfromSDSS=False)
# #tools.plot_JPLUS_results.plot_MEDspectra(new_jpl, mask=~qsomask)
# sys.exit()



#----------------------  FWHM HISTOGRAM  ----------------------#
# plt.hist(new_jpl['fwhm'], range=(0,10), bins=20, histtype='step')
# plt.xlabel('FWHM  [pixels]')
# plt.savefig(setup['plots']+'FWHM_histogram.png')
# plt.show()
# plt.close()
# sys.exit()



#----------------------  Lya LUMINOSITIES  ----------------------#
fluxlya, fcont, flux_nb, fcont_error = tools.threeFM.jplus3FM(jpl)
comov_table = setup['tools'] + 'comoving_distance.txt'
z, eta = np.genfromtxt(comov_table, skip_header=1, unpack=True)
line = interp1d(z, eta, kind='nearest', bounds_error=None, fill_value='extrapolate')
dl = line(setup['zmean'])*(1+setup['zmean'])                   #now in Mpc, all at same redshift
#dl = line(jpl['redshift'][:])*(1+jpl['redshift'][:])   #now in Mpc, redshift from leph.
dl = dl*(1.e6)*3.0856e18              #now in cm 
lumilya = 4.*np.pi*fluxlya*(dl**2.)   #now in erg/s

jpl['lya_lumin'] = lumilya



#-----------------------  DATA SAVING  -----------------------#
#sys.exit()

all_final = jplus.tools.select_object(jpl, candidates)
final = jplus.tools.select_object(jpl, over_SN)
dd.io.save(setup['final_catalog'], all_final)
dd.io.save(setup['final_catalog']+'_overSN', final)
dd.io.save(setup['flagged_catalog'], jpl)

matrix = np.empty( (8, len(jpl['coords'][candidates,0])) )
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
np.savetxt(setup['results']+'final_list_22-05-2017', matrix.T, fmt=fmt_string, header=head)


matrix = np.empty( (2, len(jpl['coords'][candidates,0])) )
matrix[0,:] = jpl['coords'][candidates,0]
matrix[1,:] = jpl['coords'][candidates,1]
head = 'ra        dec'
fmt_string = '%3.5f %3.5f'
np.savetxt(setup['results']+'final_list_ra-dec_22-05-2017', matrix.T, fmt=fmt_string, header=head)


matrix = np.empty( (2, len(jpl['coords'][over_SN,0])) )
matrix[0,:] = jpl['coords'][over_SN,0]
matrix[1,:] = jpl['coords'][over_SN,1]
head = 'ra        dec'
fmt_string = '%3.5f %3.5f'
np.savetxt(setup['results']+'final_list_ra-dec_overSN_22-05-2017', matrix.T, fmt=fmt_string, header=head)
sys.exit()


#--------------- PLOT LUMINOSITY FUNCTION ---------------#
n_tiles = 53
jplus_area = (1.4**2.)*n_tiles #now in sq.degrees
sphere_area = 4*np.pi/((np.pi/180.)**2.) #now in sq.degrees
fsky = jplus_area/sphere_area

propVol = False
tools.plot_JPLUS_results.plot_lumiFunc(fsky, line, propVol)

sys.exit()
