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
    uPrint('loading sdss galaxies... ', appendix=' ')
    sdss_gal = dd.io.load(setup['sdss_gal'])  #load sdss galaxies catalogue
    uPrint('data loading....      ', appendix='done')
    
if setup['load_sdssQSO'] == True:
    uPrint('loading sdss quasars... ', appendix=' ')
    sdss_qso = dd.io.load(setup['sdss_qso'])  #load sdss quasars catalogue
    qso_zmask = ((sdss_qso['redshift'] > setup['zmin']) & (sdss_qso['redshift'] < setup['zmax']))
    sdss_qso_z = jplus.tools.select_object(sdss_qso, qso_zmask)
    uPrint('data loading....      ', appendix='done')
    
if setup['load_gaia'] == True:
    uPrint('loading gaia... ', appendix=' ')
    gaia = dd.io.load(setup['gaia'])          #load gaia stars catalogue
    uPrint('data loading....      ', appendix='done')
    
if setup['morpho_plot'] == True or setup['morpHisto_plot'] == True:
    uPrint('loading tile_image info... ', appendix=' ')
    images = dd.io.load(setup['img_input'])    #load tile_image infos catalogue
    uPrint('data loading....           ', appendix='done')



#######################################################
#                                                     #
#                 WORK ON MOCKS DATA                  #
#                                                     #
#######################################################
if setup['load_mock_in'] == True:
    from mock_analysis import *    
    uPrint('loading mock data... ', appendix=' ')
    all_mock = tools.data_reader.mock_reader(setup['mock_input'], lines=True, emitters=True)
    uPrint('data loading....      ', appendix='done')
    uPrint('working on mock data... ', appendix=' ')
    mock = analyse_jplus_mock(all_mock)
    uPrint('working on mock data    ', appendix='done')
else:
    mock = []



#######################################################
#
#                WORK ON JPLUS DATA
#
#######################################################

uPrint('working on jplus data... ', appendix=' ')

#-------------  INITIAL MAG MASK  ------------#
mag_mask =  ((first_jpl[filt[0]][:,0] > 0) & \
             (first_jpl[filt[1]][:,0] > 0) & \
             (first_jpl[filt[2]][:,0] > 0) & \
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
extdness = mumax - jpl['rJAVA'][:,0]
    
'''    
#-------------  CROSS-MATCHES  ------------#
#-- I need this section for tile-by-tile --#
#--    plots in the following analysis   --#
uPrint('working on jplus data: cross-matching...', appendix=' ')

# SDSS GALAXIES
if setup['load_sdssGal'] == True:
    djs, ijs = ang_xmatch(jpl['coords'], sdss_gal['coords'], max_distance=1/3600.)
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
    
# GAIA STARS
if setup['load_gaia'] == True:
    djg, ijg = ang_xmatch(jpl['coords'], gaia['coords'], max_distance=1/3600.)
    mask_jg = (djg != np.inf)
else:
    mask_jg = np.ones(len(jpl['coords'][:,0]), dtype=bool)
'''

    

#----------------------------------------------------------------------------#
#------------------------  TILE-by-TILE ANALYSIS  ---------------------------#
#----------------------------------------------------------------------------#
tile=[]; ra=[]; dec=[]; flux_cont=[]; flux_cont_error=[]; cont_mag=[]; cont_mag_error=[]
jp_morph = []; sdss_morph = []; qso_morph = []; qso_morph_z = []; gaia_morph = []; ext=[[],[]]

for i in set(jpl['tile_id']):
    uPrint('working on jplus data tile by tile... now: ',appendix=str(int(i)))
    tilemask = (jpl['tile_id'] == i)
    # if (i != 2521):
    #    continue
    # if (i == 2848) or (i == 4782) or (i == 2519) or (i == 2521):# and (i != 3314):
    #    continue

    #-----------  JPLUS SIGMA LINE  -----------#
    continuum = np.array([mag_cont[tilemask], mag_cont_error[tilemask]]).T
    narr_band = np.array([jpl[filt[0]][tilemask,0], jpl[filt[0]][tilemask,1]]).T
    line = tools.sigma_line.sigmaline(continuum, narr_band)
    sigma_mask = ( mag_color > line(jpl[filt[0]][:,0]) )
    
    #---------- NRW-BAND S/N RATIO and MAG CUTS ----------#
    mcut = tools.signalTOnoise.jplus_SNratio(jpl, band=filt[0], mask=tilemask)
    snratio_mask = (jpl[filt[0]][:,0] < mcut)

    #---------- COLOR-COLOR SELECTION: ALWAYS IN [g-NB; r-NB] COLOR SPACE ----------#
    bbc_cut = setup['cCut']
    bbuc_cut = setup['cCut']
    bbc_color = jpl['gJAVA'][:,0] - jpl[filt[0]][:,0]
    bbuc_color = jpl['rJAVA'][:,0] - jpl[filt[0]][:,0]
    colcol_mask = ((bbc_color > bbc_cut) & (bbuc_color > bbuc_cut))

    #---------- COLOR SIGNAL TO NOISE RATIO SELECTION ----------#
    if setup['col_SNcut'] == True:
        color_StoN_1 = 1./np.sqrt( (jpl[filt[1]][:,1])**2. + (jpl[filt[0]][:,1])**2. )
        color_StoN_2 = 1./np.sqrt( (jpl[filt[2]][:,1])**2. + (jpl[filt[0]][:,1])**2. )
        #color_SNmask = ((color_StoN_1 > setup['colorSN']) & (color_StoN_2 > setup['colorSN']))
        color_SNmask = ((color_StoN_2 > setup['colorSN']))
    else:
        color_SNmask = np.ones(len(jpl[filt[1]][:,0]), dtype=bool)

    #-------------- SOURCE EXTENT  --------------#
    tools.evaluate_morphology(jpl, extdness, tilemask, i)
    #sys.exit()
    
    '''
    # mask_mag_tile = ((jpl['rJAVA'][tilemask,0] > 14.) & (jpl['rJAVA'][tilemask,0] < 16.))
    
    # pxl = 0.55*0.55  # jplus pixel area in arcsec^2
    # mumax = jpl['mu_max_r'][tilemask] - 2.5*np.log10(pxl)
    # extdness = mumax - jpl['rJAVA'][tilemask,0]

    avg_ext = np.mean(extdness[mask_mag_tile])
    lim_ext = avg_ext + 0.18 #EXTREMELY UNACCURATE!!!! I think I should construct a function that
    #evaluates a linear fit on the "point-like cloud" (see morpho_plot if not clear) for each tile
    #and returns it here, where it will be compared to the values in the extdness array. In this way
    #I can evaluate the probability for a source to be compact, and based on that, a "compactness
    #mask" (I foresay, this mask will be based on a cut on the "compactness probability". Then I will
    #use this maask to select tile-independent quantities (as I did for ra and dec just below here)
    #and then use those selected tile-independent values to find my compact/extended sources after
    #the end of the loop on tiles (just as I did with ra, dec and totalmask.. see below)

    #When I implement what I wrote above I would not need anymore the following line
    ext[0].append(i); ext[1].append(lim_ext)
    

    #---------- FINAL CANDIDATES LIST FOR CURRENT TILE ----------#
    totalmask = ((tilemask) & (sigma_mask) & (snratio_mask) & (colcol_mask) & (color_SNmask))
    
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
        morplot(jpl, images, mask_gaia=mask_jg, mask_sdss=mask_js, mask_quasar=mask_jq_z,\
                mask_tile=tilemask, mask_sel=totalmask, mag_cut=mkkut, tile_num=i, MUMAX=True,\
                pntLike_line=avg_ext)
            
    if setup['cstar_plot'] == True:
        cstarPlot(jpl, mask_gaia=mask_jg, mask_sdss=mask_js, mask_quasar=mask_jq,\
                  mask_tile=tilemask, mask_sel=totalmask, mag_cut=mcut, tile_num=i)
        
    if setup['morpHisto_plot'] == True:
        psf_mask = ((images['tile_id'] == i) & (images['filter_name'] == 'rJAVA'))
        morpho_pdf = jpl['fwhm'][:]/images['psf'][psf_mask][0]
        mask_sdss = ((mask_js) & (tilemask))
        mask_quasar = ((mask_jq) & (tilemask))
        mask_quasar_z = ((mask_jq_z) & (tilemask))
        mask_gaia = ((mask_jg) & (tilemask))
        jp_morph.extend(morpho_pdf[totalmask])
        sdss_morph.extend(morpho_pdf[mask_sdss])
        qso_morph.extend(morpho_pdf[mask_quasar])
        qso_morph_z.extend(morpho_pdf[mask_quasar_z])
        gaia_morph.extend(morpho_pdf[mask_gaia])

    qso_in_g = ((mask_jg) & (mask_jq) & (tilemask)) #QSOs in Gaia

    '''
print '\n'
sys.exit()

if setup['morpHisto_plot'] == True:
    tools.plot_JPLUS_results.morpHisto(jp_morph, sdss_morph, qso_morph_z, gaia_morph)

uPrint('working on jplus data tile by tile... ', appendix='          ')
uPrint('working on jplus data tile by tile... ', appendix='done')
print 'final number of candidates: ', len(ra)



#When I implement the proper morphological selections I would not need anymore the 2 following lines
ext = np.array(ext).T
np.savetxt('tile_to_extent.txt', ext, fmt='%4d   %4.2f')


sys.exit()
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
sys.exit()
#dd.io.save(setup['jplus_candidates'], new_jpl) #save list

uPrint('formatting ad saving jplus candidates catalogue... ', appendix='done')


