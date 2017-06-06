import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


import tools.settings
setup = tools.settings.set_up()


#computes N-sigma error lines for color-magnitude diagrams
def sigmaline(cont, NB, sigma=setup['sigma'], magbin=setup['mbin'], maglim=setup['ml'], show_plt=False):
    pcent = {
        1 : 68.27,
        2 : 95.45,
        3 : 99.73,
        4 : 99.99,
        5 : 99.99994,
        6 : 99.9999998,
        7 : 99.9999999997
    }.get(sigma, ' ')
    if pcent == ' ':
        raise ValueError("\n\nPercentile value associated to %s-sigma not found"%sigma)

    mask = (cont[:,1] > 0)
    color = cont[mask,0] - NB[mask,0]
    avg = np.mean(color)
    phot_err = np.sqrt(cont[mask,1]**2. + NB[mask,1]**2.)
    
    magarr = np.arange(maglim[0], maglim[1], magbin)
    sigma_line = np.zeros(len(magarr))
    cc = 0
    for k in magarr:
        idg = ((NB[mask,0] >= k) & (NB[mask,0] < k+magbin) &\
               (color > avg-setup['width']) & (color < avg+setup['width']))
        i_err = phot_err[idg] #error of the points selected by the mask
        if len(i_err) < 3:
            sigma_line[cc] = avg
        else:
            sigma_line[cc] = np.percentile(i_err, pcent)
        cc += 1

    sigma_line = interp1d(magarr,sigma_line,kind='linear',bounds_error=None,fill_value='extrapolate')

    if show_plt == True:
        # FAST PLOT TO SHOW THE ALGORITHM PERFORMANCE
        violet = (0.6, 0.0, 0.6)
        plt.plot(NB[mask,0], color, 'og', alpha=0.3, markersize=4)
        plt.plot(magarr, sigma_line(magarr), '-', c=violet, linewidth=2)
        plt.plot( (14.,24.), (avg, avg), 'r--')
        plt.axis([13.5, 24.5, -5., 10.])
        plt.show()
        plt.close()
    
    return sigma_line

