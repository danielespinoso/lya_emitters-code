import numpy as np
import sys

import matplotlib.pyplot as plt

from tools.converters import magtoflux as MtF

import tools.settings
setup = tools.settings.set_up()


def qso_lines_dict(): #from Vanden-Berk 2001
    lines = {
        'O_VI'    : 1034.0,
        'Lya'     : 1215.67,
        'Si_IV'   : 1396.76,
        'C_IV'    : 1549.06,
        'C_III'   : 1907.7,
        'Mg_II'   : 2798.75,
        'O_II'    : 3728.48,
        'Ne_III'  : 3869.85,
        'Hd'      : 4101.0,
        'Hy'      : 4340.0,
        'Hb'      : 4862.721,
        'O_III'   : 5008.239,
        'Ha'      : 6562.8
    }
    return lines



def plot_JPLUSfilters_QSOlines():
    path = setup['data'] + 'mocks/T80Cam_JAST_TransmissionCurvesTmp_20160518/'
    head = 'T80Cam_T80_';  ext = '.tab'
    
    fname = [];  waveleng = [];  transm = []
    for i in np.arange(0,12):
        filt = tools.jplus_filter_system.jpflt(i)
        fname.append(path + head + filt + ext)
        dummy, ddummy = np.loadtxt(fname[i], skiprows=1, unpack=True)
        waveleng.append(dummy)
        transm.append(ddummy)
    waveleng = np.array(waveleng)
    transm = np.array(transm)

    '''
    reds = [0.35, 1.00, 1.45, 1.7, 2.1] #redshift where one of the qso_line_dict is in J0378 filter

    f, ax = plt.subplots(len(reds), sharex=True)
    for z, j in zip(reds, range(len(reds)) ):
        kol = (1.,0,1.) #uJAVA
        ax[j].plot(waveleng[0,:], transm[0,:], c=kol, alpha=0.5, linewidth=1.8, label='uJAVA')
        kol = (0.7,0,0.7) #J0378
        ax[j].plot(waveleng[1,:], transm[1,:], c=kol, alpha=0.5, linewidth=1.8, label='J0378')
        kol = (0.5,0,1.)  #J0395
        ax[j].plot(waveleng[2,:], transm[2,:], c=kol, alpha=0.5, linewidth=1.8, label='J0395')
        kol = (0.,0.,1.)  #J0410
        ax[j].plot(waveleng[3,:], transm[3,:], c=kol, alpha=0.5, linewidth=1.8, label='J0410')
        kol = (0.,0.8,0.8)  #J0430
        ax[j].plot(waveleng[4,:], transm[4,:], c=kol, alpha=0.5, linewidth=1.8, label='J0430')
        kol = (0.,0.6,0.6)  #gJAVA
        ax[j].plot(waveleng[5,:], transm[5,:], c=kol, alpha=0.5, linewidth=1.8, label='gJAVA')
        kol = (0.4,0.8,0.)  #J0515
        ax[j].plot(waveleng[6,:], transm[6,:], c=kol, alpha=0.5, linewidth=1.8, label='J0515')
        kol = (1.,0.5,0.)  #rJAVA
        ax[j].plot(waveleng[7,:], transm[7,:], c=kol, alpha=0.5, linewidth=1.8, label='rJAVA')
        kol = (1.,0.,0.)  #J0660
        ax[j].plot(waveleng[8,:], transm[8,:], c=kol, alpha=0.5, linewidth=1.8, label='J0660')
        kol = (0.8,0.,0.)  #iJAVA
        ax[j].plot(waveleng[9,:], transm[9,:], c=kol, alpha=0.5, linewidth=1.8, label='iJAVA')
        kol = (0.6,0.,0.)  #J0861
        ax[j].plot(waveleng[10,:], transm[10,:], c=kol, alpha=0.5, linewidth=1.8, label='J0861')
        kol = (0.3,0.,0.)  #zJAVA
        ax[j].plot(waveleng[11,:], transm[11,:], c=kol, alpha=0.5, linewidth=1.8, label='zJAVA')
        
        lin = qso_lines_dict()
        
        for i in set(lin.keys()):
            ll = lin[i]*(1+z)
            if ll > 3000. and ll < 9600.:
                ax[j].plot( (ll, ll), (0.0, 0.7), 'k--', alpha=0.8, linewidth=1.8)
                ax[j].text(ll+40.,0.65, i, color='k', fontsize=8, rotation='vertical')
                
        limits=[3000., 10000., 0.0, 0.7]
        ax[j].axis(limits)
        ax[j].text(9600.,0.45, 'z = '+str(z), color='k', fontsize=8)
    plt.show()
    plt.close()
    '''

    #--------  ONE PLOT FOR REDSHIFT  --------#
    for z in np.arange(0.0, 2.6, 0.025):
        ax = plt.subplot(111)
        kol = (1.,0,1.) #uJAVA
        ax.plot(waveleng[0,:], transm[0,:], c=kol, alpha=0.5, linewidth=1.8, label='uJAVA')
        kol = (0.7,0,0.7) #J0378
        ax.plot(waveleng[1,:], transm[1,:], c=kol, alpha=0.5, linewidth=1.8, label='J0378')
        kol = (0.5,0,1.)  #J0395
        ax.plot(waveleng[2,:], transm[2,:], c=kol, alpha=0.5, linewidth=1.8, label='J0395')
        kol = (0.,0.,1.)  #J0410
        ax.plot(waveleng[3,:], transm[3,:], c=kol, alpha=0.5, linewidth=1.8, label='J0410')
        kol = (0.,0.8,0.8)  #J0430
        ax.plot(waveleng[4,:], transm[4,:], c=kol, alpha=0.5, linewidth=1.8, label='J0430')
        kol = (0.,0.6,0.6)  #gJAVA
        ax.plot(waveleng[5,:], transm[5,:], c=kol, alpha=0.5, linewidth=1.8, label='gJAVA')
        kol = (0.4,0.8,0.)  #J0515
        ax.plot(waveleng[6,:], transm[6,:], c=kol, alpha=0.5, linewidth=1.8, label='J0515')
        kol = (1.,0.5,0.)  #rJAVA
        ax.plot(waveleng[7,:], transm[7,:], c=kol, alpha=0.5, linewidth=1.8, label='rJAVA')
        kol = (1.,0.,0.)  #J0660
        ax.plot(waveleng[8,:], transm[8,:], c=kol, alpha=0.5, linewidth=1.8, label='J0660')
        kol = (0.8,0.,0.)  #iJAVA
        ax.plot(waveleng[9,:], transm[9,:], c=kol, alpha=0.5, linewidth=1.8, label='iJAVA')
        kol = (0.6,0.,0.)  #J0861
        ax.plot(waveleng[10,:], transm[10,:], c=kol, alpha=0.5, linewidth=1.8, label='J0861')
        kol = (0.3,0.,0.)  #zJAVA
        ax.plot(waveleng[11,:], transm[11,:], c=kol, alpha=0.5, linewidth=1.8, label='zJAVA')
        
        lin = qso_lines_dict()
        
        for i in set(lin.keys()):
            ll = lin[i]*(1+z)
            if ll > 3000. and ll < 9600.:
                ax.plot( (ll, ll), (0.0, 0.7), 'k--', alpha=0.8, linewidth=1.8)
                plt.text(ll+40.,0.65, i, color='k', fontsize=10, rotation='vertical')
                
        limits=[3000., 10000., 0.0, 0.7]
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax.legend(loc='center left', fontsize = 'small', bbox_to_anchor=(1, 0.5))
        ax.axis(limits)
        plt.title('redshift: '+str(int(z*100)/100.))
        plt.show()
        plt.close()

    
    sys.exit()
