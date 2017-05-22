import numpy as np
import sys
import matplotlib.pyplot as plt

from .settings import *
setup = set_up()

sys.path.append(setup['jplus_code'])
import jplus


# pipeline's example on how to use lephare,
# adapted to my needs for the lya-emitters project
def lya_lephare_jplus(data, filt_string, name='_aaa', owtspec=True, owr_prep=False, owr_run=False):
    
    for ifilter in jplus.datasets.jplus_filter_names(only_bb=False):
        if ifilter in data:
            ind = np.isfinite(data[ifilter][:,0])
            data[ifilter][~ind,:] = [-99,-99]
            data[ifilter][data[ifilter][:,0]==99,:] = [-99,-99]
            data[ifilter][data[ifilter][:,0]==0,:] = [-99,-99]
            # print ifilter
    data['redshift'] = np.zeros(len(data['rJAVA'][:,0]))
    
    Lephare = jplus.photoz.LePhare(data, per_tile=True, outspec=owtspec, recalibration=False,\
                                   zstepcfg=[0.005,2.495,0.01], emlines=True, suffix=name,\
                                   filterflag=filt_string)

    Lephare.prepare(overwrite=owr_prep)
    jp_photoz = Lephare.run(overwrite=owr_run)
    
    #LEPHARE REDSHIFT HISTOGRAM
    st_mask = (jp_photoz['bestchi2'] != 2)
    qso_mask = (jp_photoz['bestchi2'] != 1)
    qsost_mask = ((st_mask) & (qso_mask))
    plt.hist(jp_photoz['photoz'][st_mask],range=[.01,3.0],bins=100)
    plt.axvspan(2.1, 2.31, color='r', alpha=0.3, lw=0)
    plt.xlabel('z')
    plt.ylabel('# of sources')
    plt.show();plt.close()

    data['redshift'] = jp_photoz['photoz']
    return data




    
def lya_lephare_mocks(data, filt_string, name='_mocks', owtspec=True, owr_prep=False, owr_run=False):
    
    for ifilter in jplus.datasets.jplus_filter_names(only_bb=False):
        if ifilter in data:
            ind = np.isfinite(data[ifilter][:])
            data[ifilter][~ind] = [-99]
            data[ifilter][data[ifilter][:]==99] = [-99]
            data[ifilter][data[ifilter][:]==0] = [-99]
    data['redshift'] = np.zeros(len(data['rJAVA'][:,0]))

    Lephare = jplus.photoz.LePhare(data, per_tile=False, outspec=owtspec, recalibration=False, emlines=True, suffix='_mock', filterflag=filt_string)

    Lephare.prepare(overwrite=owr_prep)
    mock_photoz = Lephare.run(overwrite=owr_run)

    st_mask = (mock_photoz['bestchi2'] != 2)
    qso_mask = (mock_photoz['bestchi2'] != 1)
    qsost_mask = ((st_mask) & (qso_mask))
    #---------- REDSHIFT HISTOGRAM ----------#
    plt.hist(mock_photoz['photoz'][:],range=[.01,3.0],bins=1000)
    plt.axvspan(2.1, 2.31, color='r', alpha=0.3, lw=0)
    plt.xlabel('z')
    plt.ylabel('# of sources')
    plt.title('Mock data')
    # plt.show();plt.close()
    
    #---------- REDSHIFT PDFs ----------#
    # plt.plot(mock_photoz['PDFZSTEP'],mock_photoz['PDF'][120,:])
    # plt.show(); plt.close()
