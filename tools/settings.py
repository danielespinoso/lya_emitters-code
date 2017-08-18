
import os
import numpy as np
from . import zl

def set_up():
    setup = {}
    #######################################################
    #                                                     #
    #                     SETTINGS                        #
    #                                                     #
    #######################################################
    
    #-------------------------------#
    #----------  FOLDERS  ----------#
    #-------------------------------#
    setup['home'] = os.path.expanduser('~')+'/'            # current home folder (automatic)
    setup['here'] = setup['home'] + 'works/lya_emitters/'  # folder where 'select_emitters.py' is
    setup['jplus_code'] = setup['home'] + 'j-plus/'        # folder where Raul's jplus pipeline is
    
    setup['plots'] = setup['here'] + 'results/plots/' # where plots will be saved
    setup['data'] = setup['home'] + 'datasets/'       # where all datasets will be retrived from
    setup['tools'] = setup['here'] + 'code/tools/'    # where all tools and useful routines are
    setup['results'] = setup['here'] + 'results/'     # where final candidates-lists will be saved
    setup['spectra'] = setup['data'] + 'sdss_spec/'   # where sdss fits spectra files are stored
    
    #setup['BPZ'] = setup['home'] + 'BPZ/'        # 'root' folder of BPZ code (no more used)
    #setup['bpz'] = setup['here'] + 'bpz_files/'  # where all BPZ plots/results will be saved (unused)
    
    
    #----------------------------------#
    #----------  PARAMETERS  ----------#
    #----------------------------------#
    setup['filters'] = ['J0395', 'gJAVA', 'rJAVA']         # filters in this order: NB, BBc, BBuc
    za, zb = zl(band=setup['filters'][0], wlength=1215.67) # redshift limits corresponding to NB filter
    setup['zmin'] = za
    setup['zmax'] = zb
    setup['zmean'] = (setup['zmin'] + setup['zmax'])/2.    # mean of zmin and zmax
    
    setup['leph_all'] = [1,1,1,1,1,1,1,1,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0]  # leph. filter string (all)
    setup['leph_bb'] = [0,0,0,0,0,1,0,1,0,1,0,1, 0,0,0,0,0, 0,0,0,0,0]   # leph. filter string (bb)
    
    setup['ml'] = [14., 24.]    # general mag limits for cuts (SQL query, color masking, SNratio...)
    setup['mbin'] = 0.25        # magnitude bin for NB-excess confidence computation
    setup['SNmbin'] = 0.2       # magnitude bin for SNratio computation
    setup['marr'] = np.arange(setup['ml'][0], setup['ml'][1], 0.25)  # array of magnitudes
    setup['sigma'] = 3.         # sigmas for NB-excess confidence (see sigma_line.py routine)
    setup['SN'] = 3.            # Signal-to-Noise to determine the NB magnitude cut
    setup['bb_SN'] = 3.         # S/N cut on BB filters in workOn_emitters.py
    setup['line_ex_sigma'] = 2. # sigmas for line-excess in other NBands (see workOn_Emitters.py)
    if setup['filters'][0] == 'J0395':
        setup['width'] = 1.5    # width for NB-excess confidence (see sigma_line.py routine)
        setup['cCut'] = 0.8     # color cut for color-color diagram selection
    if setup['filters'][0] == 'J0378':
        setup['width'] = 1.     # width for NB-excess confidence (see sigma_line.py routine)
        setup['cCut'] = 0.5     # color cut for color-color diagram selection
    setup['Ha_ex'] = 0.5        # H-alpha excess cut in mags (exclude sources with rJAVA-J660 > cut)
    setup['morph_fact'] = 2.    # factor to multiply the error-fit and set the morphological limit
    setup['tile_bord'] = 100.   # number of pixels that define the tile border (if inside the border a source gets flagged)

        
    setup['mag_type'] = 'aper'      # set 'auto' or 'aper' to select auto-mag or aper-mag(3") catalogue
    setup['jpl_overwr'] = False     # if 'True' the jplus catalogue gets downloaded again from jplus website (my catalogue in lya_emitters/datasets/ gets also overwritten). Set to "True" after UPAD updates
    setup['my_overwr'] = False      # if 'True' the jplus input catalogue in lya_emitters/datasets/ gets overwritten
    setup['method'] = '2FM'         # method for NB-excess computation (choose '3FM' otherwise)
    setup['sdssPhot'] = False       # if 'True' data_reader.py (called from select_emitters.py) x-matches jplus with sdss and substitutes sdss photometry to jplus' one
    setup['galexmask'] = True       # if 'True' only sources NOT in galex will be finally selecte
    setup['morph_sel'] = 'comp'     # set 'extd' or 'comp' to respectively select extended or compact sources during the final selection (workOn_emitters.py)
    setup['data_rels'] = 'T2'      # controls which JPLUS data release is used: 'T1' for ordinary DATA-ACCESS upad catalogue, 'EDR' for Early Data Release (subset of T1) or 'T2' for test-2 data
    setup['save_firstSel'] = True  # if 'True' saves the "first_selection" catalogue (output of select_emitters.py) in the "lya_emitters/datasets" folder
        

    #-------------------------------------------------#
    #---------- JPLUS INPUT DATASETS NAMES  ----------#
    #-------------------------------------------------#    
    # JPLUS INPUTS - first jplus catalogue as formatted by datasets/data_reader.py (input of select_emitters.py)
    if setup['sdssPhot'] == True:
        setup['jplus_input'] = setup['data'] + 'jplus/jplus_allTOr24_autoMags_upad_dual_sdssPhot_sdssMatched.h5'
    else:
        setup['jplus_input'] = setup['data'] + 'jplus/jplus_allTOr24_'+setup['mag_type']+'Mags_upad_dual'+'_'+setup['data_rels']+'.h5'
    
    # JPLUS CANDIDATES - binary files produced by select_emitters.py on which 'workON_emitters.py' operates
    if setup['mag_type'] == 'aper':
        setup['jplus_candidates'] = setup['data']+'lyman_alpha/firstSelection_'+setup['filters'][0]+'_'+setup['mag_type']+'Mag_'+setup['method']+'_'+setup['data_rels']+'_SN'+str(int(setup['SN']))+'.h5'
    if setup['data_rels'] == '':
        setup['jplus_candidates'] = setup['data']+'lyman_alpha/firstSelection_'+setup['filters'][0]+'_'+setup['mag_type']+'Mag_'+setup['method']+setup['data_rels']+'_SN'+str(int(setup['SN']))+'.h5'

    

    #---------------------------------------------------------------#
    #-----------  CHOOSE WHETHER TO LOAD OTHER DATASETS  -----------#
    #----------  TO CROSSMATCH THEM IN workOn_emitters.py  ---------#
    #---------------------------------------------------------------#
    
    setup['load_mock_in'] = False   # if 'True' mock complete catalog gets loaded (in mock_analysis and select_emitters.py)
    setup['load_gaia'] = False      # if 'True' gaia stars get loaded
    setup['load_lqac'] = True       # if 'True' LQAC quasars get loaded
    setup['load_rafa'] = False      # if 'True' Rafa's HII regions get loaded
    
    setup['load_sdssGal'] = True    # if 'True' sdss galaxies catalogue gets loaded
    setup['load_sdssQSO'] = True    # if 'True' sdss quasars catalogue gets loaded
    setup['load_sdssSTAR'] = True   # if 'True' sdss stars catalogue gets loaded

    setup['load_mock_cd'] = False   # if 'True' mock candidates get loaded

        
    #-------------------------------------------#
    #---------- OTHER DATASETS NAMES  ----------#
    #-------------------------------------------#

    # MOCK INPUT - first mock catalogue as formatted by datasets/data_reader.py
    setup['mock_input'] = setup['data']+'mocks/jplus_mocks_lines_allfilt_'+ setup['filters'][0]+'emitters_linexcess.h5'
    setup['mock_candidates'] = setup['data'] + 'mocks/candidates_mocks_'+setup['filters'][0]+'_zSelected.h5'
    
    # GAIA INPUT - only in jplus footprint
    setup['gaia'] = setup['data'] + 'gaia/gaia_jpFoot.h5'
    
    # LQAC INPUT - Large Astrometric Quasar Catalogue
    setup['lqac'] = setup['data'] + 'lqac/lqac3.h5'
    
    # RAFA INPUT - catalogue of HII regions in all jplus data (by Rafael Logronho)
    setup['rafa'] = setup['data'] + 'jplus/rafa_HII.h5'
    
    # SDSS INPUTS - as downloaded from CasJobs
    setup['sdss_allphot'] = setup['data'] + 'sdss/sdss_photgalaxies_aperMags.h5'  # all sdss photometry
    setup['sdss_gal'] = setup['data'] + 'sdss/sdss_galaxies_modelMags.h5'         # sdss galaxies
    setup['sdss_qso'] = setup['data'] + 'sdss/sdss_qso_modelMags.h5'              # sdss quasars
    setup['sdss_str'] = setup['data'] + 'sdss/sdss_stars_modelMags.h5'            # sdss stars

    # GALEX INPUT - as downloaded from CasJobs (I selected magnitude fields a bit "by chance")
    setup['galex_in'] = setup['data']+'galex/galex_photoObjs_autoMags.h5'

    # TILE_IMAGE INFO INPUT
    setup['img_input'] = setup['data']+'jplus/jplus_tile_image_(03-05-2017).h5'


        
    #--------------------------------------#
    #---------- OUTPUT DATASETS  ----------#
    #--------------------------------------#
    fold = setup['results']+setup['filters'][0]+'_candidates/'+setup['data_rels']+'_SN'+str(int(setup['SN']))+'/'+setup['morph_sel']+'/'
    if not os.path.exists(fold[:-5]): os.makedirs(fold[:-5])
    if not os.path.exists(fold): os.makedirs(fold)
    
    # FINAL CANDIDATES OUTPUTS - ascii tables of selected candidates' parameters
    setup['flag_catalog'] =  fold[:-5]+'candidates_'+setup['method']+'_'+setup['filters'][0]+'_'+setup['mag_type']+'Mags_'+'FLAGS'+'_'+setup['data_rels']+'.h5'
    setup['final_catalog'] = fold+'candidates_'+setup['method']+'_'+setup['filters'][0]+'_'+setup['mag_type']+'Mags_'+setup['morph_sel']+'_'+setup['data_rels']+'.h5'
    setup['final_list'] =    fold+'candidates_'+setup['method']+'_'+setup['filters'][0]+'_'+setup['mag_type']+'Mags_'+setup['morph_sel']+'_'+setup['data_rels']+'_allpars'+'.txt'
    setup['final_radec'] =   fold+'candidates_'+setup['method']+'_'+setup['filters'][0]+'_'+setup['mag_type']+'Mags_'+setup['morph_sel']+'_'+setup['data_rels']+'_RaDec'+'.txt'
    setup['radec_overSN'] =  fold+'candidates_'+setup['method']+'_'+setup['filters'][0]+'_'+setup['mag_type']+'Mags_'+setup['morph_sel']+'_'+setup['data_rels']+'_RaDec_overSN'+'.txt'
    setup['guillaume'] =     fold[:-5]+'Guillaume_'+setup['filters'][0]+'candidates_'+setup['data_rels']+'data_'+setup['mag_type']+'Mags_'+'FLAGS'+'.h5'
        
    # GALEX OUTPUT - list of jplus-galex x-matched sources. 'workON_emitters.py' will exclude them.
    setup['galex_list'] = setup['jplus_candidates'][:(len(setup['jplus_candidates'])-3)]+'_galexCROSSMATCHEDsources.csv'

    

    #---------------------------------------#
    #-----------  PLOTS CONTROL  -----------#
    #---------------------------------------#
    
    setup['tile_plot'] = False      # if 'True' activate tile by tile plots (see select_emitters.py)
    setup['morpho_plot'] = False    # if 'True' activate morphology plots (see select_emitters.py)
    setup['cstar_plot'] = False     # if 'True' activate CLASS_STAR plots (see select_emitters.py)

    setup['plot_mock'] = False      # if 'True' mock galaxies are plotted in color-color plot --> select_emitters.py
    setup['plot_gaia'] = True      # if 'True' gaia stars are plotted in 'tile_plots' (color-magnitude and color-color) --> select_emitters.py
    setup['plot_sdssGal'] = False   # if 'True' sdss galaxies are plotted in 'tile_plots' (color-magnitude and color-color) --> select_emitters.py
    setup['plot_sdssQSO'] = False   # if 'True' sdss quasars are plotted in 'tile_plots' (color-magnitude and color-color) --> select_emitters.py
    setup['plot_sdssSTAR'] = False  # if 'True' sdss stars are plotted in 'tile_plots' (color-magnitude and color-color) --> select_emitters.py
        
    return setup
