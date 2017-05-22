import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import tools.settings
setup = tools.settings.set_up()


# calculator of the mean S/N ratio of a selection
# of jplus data (given a filter in which to operate)
def jplus_SNratio(data, mask, band='J0395', SNcut=setup['SN'], magbin=setup['SNmbin'],\
                  maglim=setup['ml']):
    sn_nb = []
    xmag = []
    for k in np.arange(maglim[0], maglim[1], magbin):
        if len(mask) == 0:
            cut = ((data[band][:,0] >= k) & (data[band][:,0] < k+magbin))
        else:
            cut = ((mask) & (data[band][:,0] >= k) & (data[band][:,0] < k+magbin))
        if (len(data[band][cut,1]) < 3):
            dummy = np.nan
        else:
            dummy = np.mean(1./data[band][cut,1]) #average S/N in the magnitude bin
        sn_nb.append(dummy)
        centre = k+(magbin/2.)
        xmag.append(centre)
    xmag = np.array(xmag)
    sn_nb = np.array(sn_nb)
    tt =  np.isnan(sn_nb)

    # when there are too little points the mag_vs_SNratio relation (blue thick line in plot) becomes
    # too noisy and mag_cut is difficult to compute. To avoid this, I perform two consecutive
    # interpolations on mag_vs_SNratio relation (2nd one is cyan thin line). Plot to understand.
    g = interp1d(xmag[~tt], sn_nb[~tt], kind='linear')
    magsNB = np.arange(min(xmag[~tt]), max(xmag[~tt]), 1.9*magbin)
    snsNB = g(magsNB)
    f = interp1d(snsNB, magsNB, kind='linear', bounds_error=None, fill_value='extrapolate')
    mag_cut = f(SNcut)

    # plt.plot(data[band][mask,0], (1./data[band][mask,1]), 'or', markersize=1)
    # plt.plot(xmag, sn_nb, 'b-', linewidth=4)
    # plt.plot(magsNB, snsNB, 'r-', linewidth=1)
    # ssNB = np.arange(0, 25, 0.05)
    # mmNB = f(ssNB)
    # plt.plot(mmNB, ssNB, 'c-', linewidth=1)
    # plt.plot( (16.,26.), (SNcut, SNcut), 'g-', linewidth=1)
    # plt.axis([16., 26., -0.5, 15.])
    # plt.show()
    # plt.close()

    return mag_cut
