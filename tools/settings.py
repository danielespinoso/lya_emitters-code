

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
    setup['data'] = setup['here'] + 'datasets/'       # where all datasets will be retrived from
    setup['tools'] = setup['here'] + 'code/tools/'    # where all tools and useful routines are
    setup['results'] = setup['here'] + 'results/'     # where final candidates-lists will be saved
    
    #setup['BPZ'] = setup['home'] + 'BPZ/'        # 'root' folder of BPZ code (no more used)
    #setup['bpz'] = setup['here'] + 'bpz_files/'  # where all BPZ plots/results will be saved (unused)
    
    
    #----------------------------------#
    #----------  PARAMETERS  ----------#
    #----------------------------------#
    setup['filters'] = ['J0378', 'uJAVA', 'gJAVA']         # filters in this order: NB, BBc, BBuc
    za, zb = zl(band=setup['filters'][0], wlength=1215.67) # redshift limits corresponding to NB filter
    setup['zmin'] = za
    setup['zmax'] = zb
    setup['zmean'] = (setup['zmin'] + setup['zmax'])/2.    # mean of zmin and zmax
    
    setup['leph_all'] = [1,1,1,1,1,1,1,1,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0]  # leph. filter string (all)
    setup['leph_bb'] = [0,0,0,0,0,1,0,1,0,1,0,1, 0,0,0,0,0, 0,0,0,0,0]   # leph. filter string (bb)
    
    setup['ml'] = [14., 24.]    # general mag limits for cuts (SQL query, color masking, SNratio...)
    setup['mbin'] = 0.25        # magnitude bin for NB-excess confidence computation
    setup['SNmbin'] = 0.2      # magnitude bin for SNratio computation
    setup['marr'] = np.arange(setup['ml'][0], setup['ml'][1], 0.25)  # array of magnitudes
    setup['sigma'] = 3.         # sigmas for NB-excess confidence (see sigma_line.py routine)
    setup['SN'] = 3.            # Signal-to-Noise to determine the NB magnitude cut
    setup['colorSN'] = 3.       # Signal-to-Noise on colors (initial selection in select_emitters.py)
    setup['line_ex_sigma'] = 2. # sigmas for line-excess in other NBands (see workOn_Emitters.py)
    if setup['filters'][0] == 'J0395':
        setup['width'] = 1.5    # width for NB-excess confidence (see sigma_line.py routine)
        setup['cCut'] = 0.8     # color cut for color-color diagram selection
    if setup['filters'][0] == 'J0378':
        setup['width'] = 1.     # width for NB-excess confidence (see sigma_line.py routine)
        setup['cCut'] = 0.5     # color cut for color-color diagram selection
    setup['Ha_ex'] = 0.5        # H-alpha excess cut in mags (exclude sources with rJAVA-J660 > cut)

        
    setup['mag_type'] = 'auto'  # set 'auto' or 'aper' to select auto-mag or aper-mag(3") catalogue
    setup['method'] = '2FM'     # method for NB-excess computation (choose '3FM' otherwise)
    setup['sdssPhot'] = False   # if 'True' x-matches jplus-sdss + uses sdss photometry inst. of jplus
    setup['zmask'] = False      # if 'True' only zmin < z < zmax sources will be finally selected
    setup['galexmask'] = False  # if 'True' only sources NOT in galex will be finally selected
    setup['col_SNcut'] = False  # if 'True' activates the color-SN cut on data (select_emitters.py)

        
    #-------------------------------------------#
    #---------- INPUT DATASETS NAMES  ----------#
    #-------------------------------------------#
    
    # JPLUS INPUTS - first jplus catalogue as formatted by datasets/data_reader.py
    if setup['sdssPhot'] == True:
        setup['jplus_input'] = setup['data'] +\
                               'jplus_alltoR24_autoMags_upad_dual_sdssPhot_sdssMatched.h5'
    else:
        if setup['mag_type'] == 'aper':
            setup['jplus_input'] = setup['data'] + 'jplus_alltoR24_aperMags_upad_dual.h5'
        elif setup['mag_type'] == 'auto':
            setup['jplus_input'] = setup['data'] + 'jplus_alltoR24_autoMags_upad_dual_jplusPhot.h5'
    
    # JPLUS CANDIDATES - binary files on which 'workON_emitters.py' operates
    if setup['mag_type'] == 'aper':
        setup['jplus_candidates'] = setup['data']+'candidates_jplus_aperMag_'+setup['filters'][0]+\
                                    '_'+setup['method']+'.h5'
    elif setup['mag_type'] == 'auto':
        setup['jplus_candidates'] = setup['data']+'candidates_jplus_autoMag_'+setup['filters'][0]+\
                                    '_'+setup['method']+'.h5'
    
    # SDSS INPUTS - as downloaded from CasJobs
    setup['sdss_allphot'] = setup['data'] + 'sdss_photgalaxies_aperMags.h5'  # all sdss photometry
    setup['sdss_gal'] = setup['data'] + 'sdss_galaxies_aperMags.h5'          # all sdss galaxies
    setup['sdss_qso'] = setup['home'] + 'jplus_data/sdss_qso_aperMags.h5'    # all sdss quasars
    
    # LQAC INPUT - Large Astrometric Quasar Catalogue
    setup['lqac'] = setup['data'] + 'lqac/lqac3.h5'
    
    # GAIA INPUT - only in jplus footprint
    setup['gaia'] = setup['data'] + 'gaia_jpFoot.h5'

    # GALEX INPUT - as downloaded from CasJobs (I selected magnitude fields a bit "by chance")
    setup['galex_in'] = setup['data']+'galex_photoObjs_autoMags.h5'

    # MOCK INPUT - first mock catalogue as formatted by datasets/data_reader.py
    setup['mock_input'] = setup['data']+'jplus_mocks_lines_allfilt_'+\
                          setup['filters'][0]+'emitters_linexcess.h5'
    setup['mock_candidates'] = setup['data'] + 'candidates_mocks_'+setup['filters'][0]+'_zSelected.h5'

    # TILE_IMAGE INFO INPUT
    setup['img_input'] = setup['data']+'jplus_tile_image_(03-05-2017).h5'
    
    # RAFA INPUT - catalogue of HII regions in all jplus data (by Rafael Logronho)
    setup['rafa'] = setup['data'] + 'rafa_HII.h5'

    

    #------------------------------------------------#
    #-----------  LOAD ADDITIONAL INPUTS  -----------#
    #------------------------------------------------#
    
    setup['load_sdssGal'] = True  # if 'True' sdss galaxies catalogue gets loaded (everywhere)
    setup['load_sdssQSO'] = True  # if 'True' sdss quasars catalogue gets loaded (everywhere)
    setup['load_mock_in'] = False  # if 'True' mock complete catalog gets loaded (select_emitters.py)
    setup['load_mock_cd'] = False  # if 'True' mock candidates get loaded (workOn_emitters.py)
    setup['load_gaia'] = False     # if 'True' gaia stars get loaded (for morphology analysis)
    setup['load_lqac'] = True     # if 'True' gaia stars get loaded (for morphology analysis)
    setup['load_rafa'] = True     # if 'True' gaia stars get loaded (for morphology analysis)
        

        
    #--------------------------------------#
    #---------- OUTPUT DATASETS  ----------#
    #--------------------------------------#
    
    # FINAL CANDIDATES OUTPUTS - ascii tables of selected candidates' parameters
    if setup['zmask'] == True:
        if setup['mag_type'] == 'auto':
            setup['final_catalog'] = setup['results'] + 'candidates_catalog_autoMags_'+\
                                     setup['filters'][0]+'_'+setup['method']+'_zSelected.h5'
            setup['flagged_catalog'] = setup['results'] + 'candidates_FLAGcatalog_autoMags_'+\
                                       setup['filters'][0]+'_'+setup['method']+'_zSelected.h5'
            setup['final_list'] = setup['results'] + 'candidates_allparams_autoMags_'+\
                                  setup['filters'][0]+'_'+setup['method']+'_zSelected.txt'
            setup['final_radec'] = setup['results'] + 'candidates_RaDec_autoMags_'+\
                                   setup['filters'][0]+'_'+setup['method']+'_zSelected.txt'
        elif setup['mag_type'] == 'aper':
            setup['final_catalog'] = setup['results'] + 'candidates_catalog_aperMags_'+\
                                     setup['filters'][0]+'_'+setup['method']+'_zSelected.h5'
            setup['flagged_catalog'] = setup['results'] + 'candidates_FLAGcatalog_aperMags_'+\
                                       setup['filters'][0]+'_'+setup['method']+'_zSelected.h5'
            setup['final_list'] = setup['results'] + 'candidates_allparams_aperMags_'+\
                                  setup['filters'][0]+'_'+setup['method']+'_zSelected.txt'
            setup['final_radec'] = setup['results'] + 'candidates_RaDec_aperMags_'+\
                                   setup['filters'][0]+'_'+setup['method']+'_zSelected.txt'
    else:
        if setup['mag_type'] == 'auto':
            setup['final_catalog'] = setup['results'] + 'candidates_catalog_autoMags_'+\
                                     setup['filters'][0]+'_'+setup['method']+'.h5'
            setup['flagged_catalog'] = setup['results'] + 'candidates_FLAGcatalog_autoMags_'+\
                                       setup['filters'][0]+'_'+setup['method']+'.h5'
            setup['final_list'] = setup['results'] + 'candidates_allparams_autoMags_'+\
                                  setup['filters'][0]+'_'+setup['method']+'.txt'
            setup['final_radec'] = setup['results'] + 'candidates_RaDec_autoMags_'+\
                                   setup['filters'][0]+'_'+setup['method']+'.txt'
        elif setup['mag_type'] == 'aper':
            setup['final_catalog'] = setup['results'] + 'candidates_catalog_aperMags_'+\
                                     setup['filters'][0]+'_'+setup['method']+'.h5'
            setup['flagged_catalog'] = setup['results'] + 'candidates_FLAGcatalog_aperMags_'+\
                                       setup['filters'][0]+'_'+setup['method']+'.h5'
            setup['final_list'] = setup['results'] + 'candidates_allparams_aperMags_'+\
                                  setup['filters'][0]+'_'+setup['method']+'.txt'
            setup['final_radec'] = setup['results'] + 'candidates_RaDec_aperMags_'+\
                                   setup['filters'][0]+'_'+setup['method']+'.txt'
            
    # GALEX OUTPUT - list of jplus-galex x-matched sources. 'workON_emitters.py' will exclude them.
    setup['galex_list'] = setup['jplus_candidates'][:(len(setup['jplus_candidates'])-3)] +\
                          '_galexCROSSMATCHEDsources.csv'



    #---------------------------------------#
    #-----------  PLOTS CONTROL  -----------#
    #---------------------------------------#
    
    setup['tile_plot'] = False        # if 'True' activate tile by tile plots (see select_emitters.py)
    setup['morpho_plot'] = False      # if 'True' activate morphology plots (see select_emitters.py)
    setup['cstar_plot'] = False       # if 'True' activate CLASS_STAR plots (see select_emitters.py)
    setup['morpHisto_plot'] = False   # if 'True' activate morph. histogram (see select_emitters.py)
    
    setup['plot_sdssGal'] = True      # if 'True' sdss galaxies are plotted in plots (col-col & morph)
    if setup['plot_sdssGal'] == True:
        setup['load_sdssGal'] = True  # forces to load the dataset
    
    setup['plot_sdssQSO'] = True      # if 'True' sdss quasars are plotted in plots (col-col & morph)
    if setup['plot_sdssQSO'] == True:
        setup['load_sdssQSO'] = True  # forces to load the dataset
    
    setup['plot_mock'] = False        # if 'True' mock galaxies are plotted in color-color plot
    if setup['plot_mock'] == True:
        setup['load_mock'] = True     # forces to load the dataset
    
    setup['plot_gaia'] = True        # if 'True' gaia stars are plotted in morphology plot
    if setup['plot_gaia'] == True:
        setup['load_gaia'] = True     # forces to load the dataset
        
    return setup
