import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import tools.settings
setup = tools.settings.set_up()


# calculator of the mean S/N ratio of a selection
# of jplus data (given a filter in which to operate)
def jplus_SNratio(data, mask, band='J0395', SNcut=setup['SN'], magbin=setup['SNmbin'], maglim=setup['ml'], show_plt=False):
    sn_nb = []
    xmag = []
    for k in np.arange(maglim[0], maglim[1], magbin):
        if len(mask) == 0:
            cut = ((data[band][:,0] >= k) & (data[band][:,0] < k+magbin))
            mask = np.ones(len(data['coords'][:,0]), dtype=bool)
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

    f = interp1d(sn_nb[~tt], xmag[~tt], kind='linear', bounds_error=None, fill_value='extrapolate')
    mag_cut = f(SNcut)

    # Test to ensure that there are enough points to compute a good S/N interpolation
    test = ( (1./data[band][:,1]) < SNcut)
    if len(data[band][test,1]) < 20:   # 20 is arbitrary.
        print '\nS/N plot might be wrongly evaluated. Check plot and tools/signalTOnoise.py!'
        mag_cut = 24. #EXTREMELY UNACCURATE!!!
        
    if show_plt == True:
        fig = plt.figure(figsize=(10,7))
        plt.plot(data[band][mask,0], (1./data[band][mask,1]), 'or', markersize=2, alpha=0.4, label='current tile data')
        plt.plot(xmag, sn_nb, 'b-', linewidth=2.5, label='S/N average per bin')
        plt.plot( (16.,26.), (SNcut, SNcut), 'g-', linewidth=1.5, label='S/N cut = '+str(int(SNcut)))
        plt.axis([16., 26., -0.5, 15.])
        plt.title('S/N interpolation test-plot')
        plt.xlabel(band+'  [mags]')
        plt.ylabel('S/N')
        plt.legend()
        plt.show()
        plt.close()

    return mag_cut
