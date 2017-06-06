import sys
import os
import deepdish as dd
import numpy as np
import matplotlib.pylab as plt

import tools

setup = tools.set_up()
filt = setup['filters']
from tools.screen_info import update_print as uPrint
from tools.converters import fluxtomag as FtoM

if setup['tile_plot'] == True:
    from tools.plot_by_tile import plot_colorMag_bytile as CM
    from tools.plot_by_tile import plot_ColCol_bytile as CC
if setup['morpho_plot'] == True:
    from tools.plot_by_tile import plot_morpho_bytile as morplot
if setup['cstar_plot'] == True:
    from tools.plot_by_tile import plot_CSTAR_bytile as cstarPlot
    
sys.path.append(setup['jplus_code'])
import jplus
from jplus.tools import crossmatch_angular as ang_xmatch



#######################################################
#                                                     #
#                  LOADING DATASETS                   #
#   MORE DETAILS IN '-->/code/datasets/README.txt'    #
#                                                     #
#######################################################

first_jpl = tools.data_reader.jplus_reader(setup['jplus_input'])

if setup['load_sdssGal'] == True:
    sdss_gal = tools.get_sdss_galaxies()  #loads sdss galaxies catalogue
    
if setup['load_sdssQSO'] == True:
    sdss_qso = tools.get_sdss_qso()       #loads sdss quasars catalogue
    qso_zmask = ((sdss_qso['zspec'] > setup['zmin']) & (sdss_qso['zspec'] < setup['zmax']))
    sdss_qso_z = jplus.tools.select_object(sdss_qso, qso_zmask)
    
if setup['load_sdssSTAR'] == True:
    sdss_str = tools.get_sdss_stars()     #loads sdss stars catalogue
    
if setup['load_gaia'] == True:
    uPrint('loading gaia... ', appendix=' ')
    gaia = dd.io.load(setup['gaia'])          #load gaia stars catalogue
    uPrint('loading gaia....      ', appendix='done')



#######################################################
#                                                     #
#                 WORK ON MOCKS DATA                  #
#                                                     #
#######################################################
if setup['load_mock_in'] == True:
    from mock_analysis import *    
    uPrint('loading mock data... ', appendix=' ')
    all_mock = tools.data_reader.mock_reader(setup['mock_input'], lines=True, emitters=True)
    uPrint('loading mock data....      ', appendix='done')
    uPrint('working on mock data... ', appendix=' ')
    mock = analyse_jplus_mock(all_mock)
    uPrint('working on mock data    ', appendix='done')
else:
    mock = []



#######################################################
#                                                     #
#            PRELIMINARY WORK ON JPLUS DATA           #
#                                                     #
#######################################################
uPrint('working on jplus data... ', appendix=' ')

#-------------  INITIAL MAG MASK  ------------#
mag_mask =  ((first_jpl[filt[0]][:,0] > 0) & \
             (first_jpl[filt[1]][:,0] > 0) & \
             (first_jpl[filt[2]][:,0] > 0) & \
             (first_jpl['mag_auto_r'][:] > setup['ml'][0]) & \
             (first_jpl['rJAVA'][:,0] > setup['ml'][0]) & (first_jpl['rJAVA'][:,0] < setup['ml'][1]))
jpl = jplus.tools.select_object(first_jpl, mag_mask)


#--------  NARROW BAND EXCESS COMPUTATION  -------#
if (setup['method'] == '3FM'):
    fline, fcont, flux_nb, fcont_error = tools.threeFM.jplus3FM(jpl)
    mag_cont, mag_cont_error = FtoM(fcont, band=filt[0], dflux=fcont_error)
    mag_color = mag_cont - jpl[filt[0]][:,0]
elif (setup['method'] == '2FM'):
    mag_cont = jpl[filt[1]][:,0]
    mag_cont_error = jpl[filt[1]][:,1]
    mag_color = jpl[filt[1]][:,0] - jpl[filt[0]][:,0]


#-------  EXTENDED-NESS PARAMETER COMPUTATION  ------#
pxl = 0.55*0.55  # jplus pixel area in arcsec^2
mumax = jpl['mu_max_r'][:] - 2.5*np.log10(pxl)
extdness = mumax - jpl['mag_auto_r'][:]
jpl['extended'] = np.zeros(len(jpl['mu_max_r']), dtype=bool)
jpl['compact'] = np.zeros(len(jpl['mu_max_r']), dtype=bool)
    

#-------------  CROSS-MATCHES  (only for tile-by-tile plots)------------#
if (setup['tile_plot'] == True) or (setup['morpho_plot'] == True) or (setup['cstar_plot'] == True):
    uPrint('working on jplus data: cross-matching...', appendix=' ')
    
    # SDSS GALAXIES
    if setup['load_sdssGal'] == True:
        djs, ijs = ang_xmatch(jpl['coords'], sdss_gal['coords'], max_distance=3/3600.)
        mask_js = (djs != np.inf)
    else:
        mask_js = np.ones(len(jpl['coords'][:,0]), dtype=bool)
        
    # SDSS QUASARS
    if setup['load_sdssQSO'] == True:
        djq, ijq = ang_xmatch(jpl['coords'], sdss_qso['coords'], max_distance=3/3600.)
        mask_jq = (djq != np.inf)
        djq_z, ijq_z = ang_xmatch(jpl['coords'], sdss_qso_z['coords'], max_distance=3/3600.)
        mask_jq_z = (djq_z != np.inf)
    else:
        mask_jq = np.ones(len(jpl['coords'][:,0]), dtype=bool)
        mask_jq_z = np.ones(len(jpl['coords'][:,0]), dtype=bool)
    
    # SDSS STARS
    if setup['load_sdssSTAR'] == True:
        djstr, ijstr = ang_xmatch(jpl['coords'], sdss_str['coords'], max_distance=3/3600.)
        mask_jstr = (djstr != np.inf)
    else:
        mask_jstr = np.ones(len(jpl['coords'][:,0]), dtype=bool)
        
    # GAIA STARS
    if setup['load_gaia'] == True:
        djg, ijg = ang_xmatch(jpl['coords'], gaia['coords'], max_distance=3/3600.)
        mask_jg = (djg != np.inf)
    else:
        mask_jg = np.ones(len(jpl['coords'][:,0]), dtype=bool)

    uPrint('working on jplus data: cross-matching...', appendix='done')


########################################################################################
#                                                                                      #
#                         TILE-by-TILE ANALYSIS ON JPLUS DATA                          #
#                                                                                      #
########################################################################################
tile=[]; ra=[]; dec=[]; flux_cont=[]; flux_cont_error=[]; cont_mag=[]; cont_mag_error=[]
jp_morph = []; sdss_morph = []; qso_morph = []; qso_morph_z = []; gaia_morph = []; ext=[[],[]]
print '\nNumber of tiles to analyse: ', len(set(jpl['tile_id']))

for i in set(jpl['tile_id']):
    uPrint('working on jplus data tile by tile... now: ',appendix=str(int(i)) )
    tilemask = (jpl['tile_id'] == i)
    # if (i == 2672):
    #   continue

    #-----------  JPLUS SIGMA LINE  -----------#
    continuum = np.array([mag_cont[tilemask], mag_cont_error[tilemask]]).T
    narr_band = np.array([jpl[filt[0]][tilemask,0], jpl[filt[0]][tilemask,1]]).T
    line = tools.sigma_line.sigmaline(continuum, narr_band, show_plt=False)
    sigma_mask = ( mag_color > line(jpl[filt[0]][:,0]) )
    
    #---------- NRW-BAND S/N RATIO and MAG CUTS ----------#
    mcut = tools.signalTOnoise.jplus_SNratio(jpl, band=filt[0], mask=tilemask, show_plt=False)
    snratio_mask = (jpl[filt[0]][:,0] < mcut)

    #---------- COLOR-COLOR SELECTION: ALWAYS IN [g-NB; r-NB] COLOR SPACE ----------#
    bbc_cut = setup['cCut']
    bbuc_cut = setup['cCut']
    bbc_color = jpl['gJAVA'][:,0] - jpl[filt[0]][:,0]
    bbuc_color = jpl['rJAVA'][:,0] - jpl[filt[0]][:,0]
    colcol_mask = ((bbc_color > bbc_cut) & (bbuc_color > bbuc_cut))

    #-------------- SOURCE EXTENT  --------------#
    pp = tools.evaluate_morphology(jpl, extdness, tilemask, i) 
    border = setup['morph_fact'] * (pp[2] + (10.**(pp[3]*jpl['rJAVA'][:,0]+pp[4])) )\
             + (pp[0] + pp[1]*jpl['rJAVA'][:,0])
    extended_mask = ((extdness > border) & (tilemask))
    compact_mask = ((extdness <= border) & (tilemask))        
    jpl['extended'] = ((jpl['extended']) | (extended_mask))
    jpl['compact'] = ((jpl['compact']) | (compact_mask))
    
    #---------- FINAL CANDIDATES LIST FOR CURRENT TILE ----------#
    totalmask = ((tilemask) & (sigma_mask) & (snratio_mask) & (colcol_mask) )
    
    temp_tile_id = jpl['tile_id'][totalmask]
    temp_ra = jpl['coords'][totalmask,0]
    temp_dec = jpl['coords'][totalmask,1]
    tile.extend(temp_tile_id)
    ra.extend(temp_ra)
    dec.extend(temp_dec)
    
    #---------- FINAL PLOTS FOR CURRENT TILE ----------#
    if setup['tile_plot'] == True:
        # Color-Magnitude plot down here
        CM(jpl,mag_color,mask_selec=totalmask,mask_tile=tilemask,linea=line,mag_cut=mcut,tile_num=i)
        # Color-Color plot down here
        CC(jpl, mock, mask_sdss=mask_js, mask_tile=tilemask, mask_sel=totalmask, tile_num=i)
        
    if setup['morpho_plot'] == True:
        mkkut = tools.signalTOnoise.jplus_SNratio(jpl, band='rJAVA', mask=tilemask)
        morplot(jpl, mask_gaia=mask_jg, mask_sdss=mask_js, mask_quasar=mask_jq_z,\
                mask_tile=tilemask, mask_sel=totalmask, mag_cut=mkkut, tile_num=i, MUMAX=True,\
                bord_params=pp, ext_mask=extended_mask, comp_mask=compact_mask)
            
    if setup['cstar_plot'] == True:
        cstarPlot(jpl, mask_gaia=mask_jg, mask_sdss=mask_js, mask_quasar=mask_jq,\
                  mask_tile=tilemask, mask_sel=totalmask, mag_cut=mcut, tile_num=i)

    uPrint('working on jplus data tile by tile... now: ',appendix='        ' )    

#uPrint('working on jplus data tile by tile... ', appendix='          ')
uPrint('working on jplus data tile by tile... ', appendix='done')
print 'final number of candidates: ', len(ra)


#sys.exit()
#---------- FINAL CONVERSIONS ----------#
uPrint('formatting ad saving jplus candidates catalogue... ', appendix=' ')
tile = np.array(tile)
ra = np.array(ra)
dec = np.array(dec)

#---------- FORMATTING and SAVING OUTPUT ----------#
ramask = np.in1d(jpl['coords'][:,0], ra[:])
decmask = np.in1d(jpl['coords'][:,1], dec[:])

finalmask = ((ramask) & (decmask))

new_jpl = jplus.tools.select_object(jpl, finalmask)

new_jpl['redshift'] = (setup['zmean'])*np.ones(len(new_jpl['rJAVA'][:,0]))
#sys.exit()
dd.io.save(setup['jplus_candidates'], new_jpl) #save list

uPrint('formatting ad saving jplus candidates catalogue... ', appendix='done')


