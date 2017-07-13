import numpy as np
import sys
import os
import deepdish as dd

import tools
setup = tools.set_up()
import tools.screen_info as tsi

# THIS FUNCTION READS "FULL" JPLUS AND SDSS CATALOGUES (the ones found in $HOME/jplus_data/)
# AND PREPARES THEM FOR THE CODE:  $HOME/lya_emitters/select_emitters.py   (WHICH MEANS: cleaning
# invalid values, x-matching them and substituting sdss photometry to the jplus one). At the end of
# the process it saves a new jplus catalogue in $HOME/lya_emitters/datasets. This catalogue will have
# -99. instead of [zeros, nan, +99] values and will get written only if doesn't already exist)

def jplus_reader(out_name):
    if setup['jpl_overwr'] == True: # if the dataset in home/jplus_data gets overwritten, then also the catalogue in lya_emitters/dataset must be overwritten
        setup['my_overwr'] = True
    if not os.path.exists(out_name) or setup['my_overwr']==True:
        sys.path.append(setup['jplus_code'])
        import jplus
        from jplus.datasets import fetch_jplus_objects as fetch
    
        #----------------------------  (DOWN)LOADING JPLUS  ---------------------------#
        if setup['data_rels'] == 'T1' :  dab = 'upad'      #see Raul's jplus config.ini file
        if setup['data_rels'] == 'T2' :  dab = 'test2'     #see Raul's jplus config.ini file
        if setup['data_rels'] == 'EDR':  dab = 'edr'       #see Raul's jplus config.ini file
        jpl = fetch(db=dab, mag_type=setup['mag_type']+"Mags", mag_limit=[0,24], cstar_query=" ",\
                    object_name='allTOr24', overwrite=setup['jpl_overwr'], dualMode=True,\
                    filter_name="rJAVA", nchunks=20, extra_where_conds='')
        #---------------------------  CLEANING JPLUS  --------------------------------#
        for ifilter in jplus.datasets.jplus_filter_names():
            ind = np.isfinite(jpl[ifilter][:,0]) #cleaning NaNs
            jpl[ifilter][~ind,:] = [-99.0,-99.0]
            mask = (jpl[ifilter][:,0] == 0.) #cleaning zeros
            jpl[ifilter][mask,0] = -99.0
            jpl[ifilter][mask,1] = -99.0
            mask = (jpl[ifilter][:,0] == 99) #cleaning +99
            jpl[ifilter][mask,0] = -99.0
            jpl[ifilter][mask,1] = -99.0
        iind = np.isfinite(jpl['mag_auto_r'][:]) #cleaning NaNs from "mag_auto_r"
        jpl['mag_auto_r'][~iind] = -99.0
        mask = (jpl['mag_auto_r'][:] == 0.) #cleaning zeros from "mag_auto_r"
        jpl['mag_auto_r'][mask] = -99.0
        mask = (jpl['mag_auto_r'][:] == 99) #cleaning +99 from "mag_auto_r"
        jpl['mag_auto_r'][mask] = -99.0
        #----------------------------  LOADING SDSS  ---------------------------#
        if setup['sdssPhot'] == True:
            tsi.update_print('working on datasets... ', appendix='loading  sdss')
            sdss = dd.io.load(setup['sdss_allphot'])
            sdss['uJAVA'] = sdss.pop('uSDSS')
            sdss['gJAVA'] = sdss.pop('gSDSS')
            sdss['rJAVA'] = sdss.pop('rSDSS')
            sdss['iJAVA'] = sdss.pop('iSDSS')
            sdss['zJAVA'] = sdss.pop('zSDSS')
            #----------------  CROSS-MATCHING SDSS and OVERWRITING ITS BROAD BANDS  ----------------#
            tsi.update_print('working on datasets... ', appendix='x-matching catalogues')
            dist_js, ind_js = jplus.tools.crossmatch_angular(jpl['coords'], sdss['coords'],\
                                                             max_distance=1/3600.)
            mask_js = ( (dist_js != np.inf) )
            tsi.update_print('working on datasets... ', appendix='overwriting sdss photometry')
            for ifilter in jplus.datasets.jplus_filter_names(only_bb=True):
                jpl[ifilter][mask_js,0] = sdss[ifilter][ind_js[mask_js],0]
                jpl[ifilter][mask_js,1] = sdss[ifilter][ind_js[mask_js],1]
        else:
            mask_js = np.ones(len(jpl['rJAVA'][:,0]), dtype=bool)
        #------------------------  FORMATTING AND SAVING OUTPUT  -----------------------#
        new_data = {
            'filename':out_name,
            'date':jpl['date'],
            'SQL_query':jpl['SQL_query'],
            'coords':jpl['coords'][mask_js],
            'object_id':jpl['object_id'][mask_js],
            'tile_id':jpl['tile_id'][mask_js],
            'cstar':jpl['cstar'][mask_js],
            'uJAVA':jpl['uJAVA'][mask_js],
            'J0378':jpl['J0378'][mask_js],
            'J0395':jpl['J0395'][mask_js],
            'J0410':jpl['J0410'][mask_js],
            'J0430':jpl['J0430'][mask_js],
            'gJAVA':jpl['gJAVA'][mask_js],
            'J0515':jpl['J0515'][mask_js],
            'rJAVA':jpl['rJAVA'][mask_js],
            'J0660':jpl['J0660'][mask_js],
            'iJAVA':jpl['iJAVA'][mask_js],
            'J0861':jpl['J0861'][mask_js],
            'zJAVA':jpl['zJAVA'][mask_js],
            'flag_uJAVA':jpl['flag_uJAVA'][mask_js],
            'flag_J0378':jpl['flag_J0378'][mask_js],
            'flag_J0395':jpl['flag_J0395'][mask_js],
            'flag_J0410':jpl['flag_J0410'][mask_js],
            'flag_J0430':jpl['flag_J0430'][mask_js],
            'flag_gJAVA':jpl['flag_gJAVA'][mask_js],
            'flag_J0515':jpl['flag_J0515'][mask_js],
            'flag_rJAVA':jpl['flag_rJAVA'][mask_js],
            'flag_J0660':jpl['flag_J0660'][mask_js],
            'flag_iJAVA':jpl['flag_iJAVA'][mask_js],
            'flag_J0861':jpl['flag_J0861'][mask_js],
            'flag_zJAVA':jpl['flag_zJAVA'][mask_js],
            'fwhm':jpl['fwhm'][mask_js],
            'KronRad':jpl['kron_rad'][mask_js],           
            'PetroRad':jpl['petro_rad'][mask_js],
            'mag_auto_r':jpl['mag_auto_r'][mask_js],
            'mu_max_r':jpl['mu_max_r'][mask_js]
        }
        if setup['mag_type'] == 'aper':
            new_data['mag_auto_r'] = jpl['mag_auto_r'][mask_js]
            
        tsi.update_print('working on datasets...   saving  final  catalogue', appendix=' ')  
        dd.io.save(out_name, new_data)
        tsi.update_print('working on datasets...                           ', appendix='done')
    else:
        tsi.update_print('loading jplus '+setup['data_rels']+' data... ', appendix=' ')
        new_data = dd.io.load(out_name)
        tsi.update_print('loading jplus '+setup['data_rels']+' data... ', appendix='done')

    return new_data
        
        


        

#####################################################################################################
#---------------------------------------------------------------------------------------------------#
#####################################################################################################


    
#reads David's mock and writes a hdf5 dictionary output file
def mock_reader(out_name, nCone=512, csmlgy = [73.0, 0.25], filts='all', zedspace=True, lines=True,\
                emitters=True, em_filter=setup['filters'][0]):
    if not os.path.exists(out_name):
        from astropy.cosmology import FlatLambdaCDM
        from mocks.MocksCluster import readmock_chunk
        from scipy.interpolate import interp1d
        
        cosmo = FlatLambdaCDM(H0 = csmlgy[0], Om0 = csmlgy[1])
        c = 3*(10**18)  #in angstrom/s
        
        gal = np.dtype([('Type',np.int32),('Mvir',np.float32),('pos',np.float32,3),\
                        ('vel',np.float32,3),('sfr',np.float32),('BulgeMass',np.float32),\
                        ('DiskMass',np.float32),('Time',np.float32),('redshift',np.float32),
                        ('mag',np.float32),('MetalColdGas',np.float32),('ColdGas',np.float32),\
                        ('MassWeightAge',np.float32),('ObsMagDust',np.float32,12)])
        path_lines = setup['data'] + 'mocks/lines/'
        path_cont = setup['data'] + 'mocks/continuum/'
        
        if lines == True:
            appn = 'lines'
        else:
            emitters = False
            appn = 'cont'
            
        if filts == 'all':
            appn2 = '_allfilt'
        else:
            appn2 ==''
            
        ltc = 'LightCone_SA_0_'

        #------ LINE FIELDS ------#
        # Type_l = []
        # Mvir_l = []
        # pos_l = [],[],[]
        # vel_l = [],[],[]
        # sfr_l = []
        # BulMass_l = []
        # diskMass_l = []
        # Time_l = []
        # z_l = []
        # rMag_l = []
        # MetalColdGas_l = []
        # ColdGas_l = []
        # MassWeightAge_l = []
        ObsMagDust_l = [],[],[],[],[],[],[],[],[],[],[],[]
        
        #------ CONTINUUM FIELDS ------#
        # Type_c = []
        # Mvir_c = []
        # pos_c = [],[],[]
        # vel_c = [],[],[]
        # sfr_c = []
        # BulMass_c = []
        # diskMass_c = []
        # Time_c = []
        # z_c = []
        # rMag_c = []
        # MetalColdGas_c = []
        # ColdGas_c = []
        # MassWeightAge_c = []
        ObsMagDust_c = [],[],[],[],[],[],[],[],[],[],[],[]
        
        z_corr = []
        z_dist = []
        
        for i in np.arange(0,nCone,1):
            name_lines = path_lines+ltc
            a = readmock_chunk(0, 1, i, path_lines, zspace=zedspace)
            dummy = a[0]
            z_dist.extend(a[3])
            z_corr.extend(a[4])
            
            name_cont = path_cont+ltc
            b = readmock_chunk(0, 1, i, path_cont, zspace=zedspace)
            ddummy = b[0]
            
            # Type_l.extend(dummy['Type'][:]);  Type_c.extend(ddummy['Type'][:])
            # Mvir_l.extend(dummy['Mvir'][:]);  Mvir_c.extend(ddummy['Mvir'][:])
            # for k in np.arange(0,3):
            #    pos_l[k].extend(dummy['pos'][:,k]);  pos_c[k].extend(ddummy['pos'][:,k])
            #    vel_l[k].extend(dummy['vel'][:,k]);  vel_c[k].extend(ddummy['vel'][:,k])
            # sfr_l.extend(dummy['sfr'][:]);  sfr_c.extend(ddummy['sfr'][:])
            # BulMass_l.extend(dummy['BulgeMass'][:]);  BulMass_c.extend(ddummy['BulgeMass'][:])
            # diskMass_l.extend(dummy['DiskMass'][:]);  diskMass_c.extend(ddummy['DiskMass'][:])
            # Time_l.extend(dummy['Time'][:]);  Time_c.extend(ddummy['Time'][:])
            # z_l.extend(dummy['redshift'][:]);  z_c.extend(ddummy['redshift'][:])
            # rMag_l.extend(dummy['mag'][:]);  rMag_c.extend(ddummy['mag'][:])
            # MetalColdGas_l.extend(dummy['MetalColdGas'][:])
            # MetalColdGas_c.extend(ddummy['MetalColdGas'][:])
            # ColdGas_l.extend(dummy['ColdGas'][:]);  ColdGas_c.extend(ddummy['ColdGas'][:])
            # MassWeightAge_l.extend(dummy['MassWeightAge'][:])
            # MassWeightAge_c.extend(ddummy['MassWeightAge'][:])
            for k in np.arange(0,12):
                ObsMagDust_l[k].extend(dummy['ObsMagDust'][:,k])
                ObsMagDust_c[k].extend(ddummy['ObsMagDust'][:,k])
                
            tsi.update_print('reading mocks...', appendix=(i*100/(nCone-1)), percent=True)
            
        tsi.update_print('converting lists to vectors...', appendix=' ')
        # Type_l = np.array(Type_l); Mvir_l = np.array(Mvir_l)
        # pos_l = np.array(pos_l); vel_l = np.array(pos_l)
        # sfr_l = np.array(sfr_l);BulMass_l = np.array(BulMass_l);diskMass_l = np.array(diskMass_l)
        # Time_l = np.array(Time_l); z_l = np.array(z_l); rMag_l = np.array(rMag_l);
        # MetalColdGas_l = np.array(MetalColdGas_l); ColdGas_l = np.array(ColdGas_l)
        # MassWeightAge_l = np.array(MassWeightAge_l)
        ObsMagDust_l = np.array(ObsMagDust_l)
        # Type_c = np.array(Type_c); Mvir_c = np.array(Mvir_c)
        # pos_c = np.array(pos_c); vel_c = np.array(pos_c)
        # sfr_c = np.array(sfr_c);BulMass_c = np.array(BulMass_c);diskMass_c = np.array(diskMass_c)
        # Time_c = np.array(Time_c); z_c = np.array(z_c); rMag_c = np.array(rMag_c);
        # MetalColdGas_c = np.array(MetalColdGas_c); ColdGas_c = np.array(ColdGas_c)
        # MassWeightAge_c = np.array(MassWeightAge_c) 
        ObsMagDust_c = np.array(ObsMagDust_c)
        z_corr = np.array(z_corr)
        z_dist = np.array(z_dist)
        tsi.update_print('converting lists to vectors...', appendix='done')
        
        tsi.update_print('converting magnitudes from absolute to AB...', appendix=' ')
        dl = (1.+z_corr)*(z_dist/cosmo.h)
        for k in np.arange(0,12):
            ObsMagDust_l[k] = ObsMagDust_l[k] + 5*(np.log10(dl[:])+5)
            ObsMagDust_c[k] = ObsMagDust_c[k] + 5*(np.log10(dl[:])+5)
        tsi.update_print('converting magnitudes from absolute to AB...', appendix='done')
        
        #------------ SIMULATING REALISTIC MAG_ERRORS FROM JPLUS DATA ------------#
        tsi.update_print('computing mag errors...', appendix='')
        jpl = dd.io.load(setup['jplus_input'])
        jpfilters = ['uJAVA', 'J0378', 'J0395', 'J0410', 'J0430', 'gJAVA', 'J0515', 'rJAVA', 'J0660', 'iJAVA', 'J0861', 'zJAVA']
        tile_mask = (jpl['tile_id'][:] == 2525)
        mag_errs = {}
        for i, ifilt in zip(np.arange(0,12), jpfilters):
            x = jpl[ifilt][tile_mask,0]
            y = jpl[ifilt][tile_mask,1]
            line = interp1d(x, y, kind='nearest', bounds_error=None, fill_value='extrapolate')
            mag_errs[ifilt] = line(ObsMagDust_l[i])
        tsi.update_print('computing mag errors...', appendix='done')
        
        #------------- NARROW-BAND EMITTERS SELECTION MASK -------------#
        tsi.update_print('formatting data...', appendix='')
        if emitters == True:
            filt = {
                'uJAVA':0,
                'J0378':1,
                'J0395':2,
                'J0410':3,
                'J0430':4,
                'gJAVA':5,
                'J0515':6,
                'rJAVA':7,
                'J0660':8,
                'iJAVA':9,
                'J0861':10,
                'zJAVA':11,            
            }.get(em_filter, ' ')
            if filt == ' ':
                raise ValueError("\n\nFilter %s not found"%em_filter)
            excess = ObsMagDust_c[filt]-ObsMagDust_l[filt]  #continuum-line in the emission filter
            mask = (excess > 0.5)
            appn3 = '_'+em_filter+'emitters'
        else:
            mask = (ObsMagDust_c[0] > -1.e-500) #True everywhere (I hope)
            appn3 = ''
        
        #------------------------  FORMATTING AND SAVING OUTPUT  -----------------------#
        tiles = np.zeros(len(z_corr))
        tiles[:] = int(2525)
        
        archive_file = setup['data'] + 'jplus_mocks_'+appn+appn2+appn3+'_linexcess.h5'
        table = {
            'tile_id':tiles[mask],
            'NB_excess':excess[mask],
            #'type':Type_l[mask],
            #'mvir':Mvir_l[mask],
            #'coords':pos_l[:,mask].T,
            #'vel':vel_l[:,mask].T,
            #'sfr':sfr_l[mask],
            #'bulge_mass':BulMass[mask],
            #'disk_mass':diskMass[mask],
            #'time':Time[mask],
            #'redshift_uncorr':z[mask],
            #'rJAVA':rMag[mask],
            #'metallicity':MetalColdGas[mask],
            #'cold_gas':ColdGas[mask],
            #'massWeightAge':MassWeightAge[mask],
            'uJAVA':np.stack((ObsMagDust_l[0][mask],mag_errs['uJAVA'][mask]), axis=-1),
            'J0378':np.stack((ObsMagDust_l[1][mask],mag_errs['J0378'][mask]), axis=-1),
            'J0395':np.stack((ObsMagDust_l[2][mask],mag_errs['J0395'][mask]), axis=-1),
            'J0410':np.stack((ObsMagDust_l[3][mask],mag_errs['J0410'][mask]), axis=-1),
            'J0430':np.stack((ObsMagDust_l[4][mask],mag_errs['J0430'][mask]), axis=-1),
            'gJAVA':np.stack((ObsMagDust_l[5][mask],mag_errs['gJAVA'][mask]), axis=-1),
            'J0515':np.stack((ObsMagDust_l[6][mask],mag_errs['J0515'][mask]), axis=-1),
            'rJAVA':np.stack((ObsMagDust_l[7][mask],mag_errs['rJAVA'][mask]), axis=-1),
            'J0660':np.stack((ObsMagDust_l[8][mask],mag_errs['J0660'][mask]), axis=-1),
            'iJAVA':np.stack((ObsMagDust_l[9][mask],mag_errs['iJAVA'][mask]), axis=-1),
            'J0861':np.stack((ObsMagDust_l[10][mask],mag_errs['J0861'][mask]), axis=-1),
            'zJAVA':np.stack((ObsMagDust_l[11][mask],mag_errs['zJAVA'][mask]), axis=-1),
            'redshift':z_corr[mask],
            'filename':archive_file
        }
        
        table['object_id'] = np.arange(0,len(table['rJAVA']),1)+1 #object ids (from 1 to #_of_obj)
        table['coords'] = np.random.rand(len(table['rJAVA']), 2)  #random "realistic" fake-coords
        table['coords'][:,0] = 150.*table['coords'][:,0]
        table['coords'][:,1] = 35.*table['coords'][:,1]
        tsi.update_print('formatting data...', appendix='done')
        
        tsi.update_print('saving catalogue to %s...'%archive_file, appendix=' ')
        dd.io.save(archive_file, table)
        tsi.update_print('saving catalogue to %s   '%archive_file, appendix='done')
    else:
        tsi.update_print('loading mocks... ', appendix=' ')
        table = dd.io.load(out_name)
        tsi.update_print('loading mocks... ', appendix='done')

    return table





    
#cuts the un-masked part of a jplus-mock catalogue
#and returns only the masked part of it
def select_mock_object(data, position):
    data_sel = {}
    for key in data.keys():
        if key == 'sfr':
            continue
        if (type(data[key]) is not np.str_) and (type(data[key]) is not np.str):
            data_sel[key] = data[key][position]
        else:
            data_sel[key] = data[key]

    return data_sel
