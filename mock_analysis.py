import sys
import deepdish as dd
import numpy as np
import matplotlib.pylab as plt

import tools

setup = tools.set_up()
filt = setup['filters']
from tools.screen_info import update_print as uPrint
from tools.converters import fluxtomag as FtoM
from tools.lephare_lya import lya_lephare_mocks as LP_mocks

sys.path.append(setup['jplus_code'])
import jplus

def analyse_jplus_mock(data):
    #---------- REDSHIFT SELECTION ----------#
    if setup['load_mock_in'] == True:
        zmask = ((mock['redshift'][:] > setup['zmin']) & (mock['redshift'][:] < setup['zmax']))
        mock_cmask = ((zmask) & \
                      (mock[filt[1]][:,0] > 0.) & (mock[filt[2]][:,0] > 0.) &\
                      (mock[filt[0]][:,0] > 0.) & (mock[filt[0]][:,0] < 24.))
        mock = tools.data_reader.select_mock_object(mock, mock_cmask)
        #dd.io.save(setup['mock_candidates'], mock)
    else:
        mock = []
        
    '''
    #---------- LEPHARE CHECK (it should retrive the correct redshift) ----------#
    #LP_mocks(mock, setup['leph_bb'], name='_mockphz', owtspec=True, owr_prep=False, owr_run=True)
    
    #-----------  THREE FILTER METHOD ON MOCKS  ----------#
    uPrint('calculating 3FM on mocks... ', appendix=' ')
    if (setup['method'] == '3FM'):
        mfcont, mfline, mflux_nb = tools.threeFM.mocks3FM(mock)
        mmag_cont = FtoM(mfcont, band=filt[0])
        mmag_color = mmag_cont - mock[filt[0]][:,0]
        mEW = ((mflux_nb-mfcont)/mfcont)*(1.+mock['redshift']) #DOUBT: COMPUTE CONTINUUM FLUX DENSITY?
    elif (setup['method'] == '2FM'):
        mmag_cont = mock[filt[1]][:,0]
        mmag_color = mmag_cont - mock[filt[0]][:,0]
    uPrint('calculating 3FM on mocks... ', appendix='done')

    
    #-----------  PLOTTING  ----------#
    plt.hist(mmag_color, bins=60, range=(-1.,3.5), color='r', alpha=0.6, label=setup['method'])
    plt.hist(mock['NB_excess'], bins=60, range=(-1.,3.5), color='b', alpha=0.3, label='real')
    plt.xlabel(r'$\Delta$m')
    plt.legend()
    plt.show()
    plt.close()
    
    plt.plot(mock['NB_excess'], mmag_color, 'og', markersize=5, alpha=0.4)
    plt.plot( (0.5, 3.), (0.5, 3.), 'r--')
    plt.xlabel('mock real excess [mags]', fontsize=12)
    plt.ylabel('3FM excess [mags]', fontsize=12)
    plt.title(r'$\Delta$m', fontsize=12)
    #plt.axis([12., 22., 12., 22])
    plt.show()
    plt.close()

    '''
    return mock
