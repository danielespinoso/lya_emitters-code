import numpy as np
import sys

import tools.settings
from tools.converters import magtoflux as MtF
setup = tools.settings.set_up()


# Computes MEDIAN alpha and beta coefficients for a given filter
# bbc = Broad Band (filter) Contaminated (by the emission line)
# bbuc = Broad Band (filter) Uncontaminated (by the emission line)
def alphabeta(filset=setup['filters'], data_dir=setup['data']):
    path = data_dir + 'mocks/T80Cam_JAST_TransmissionCurvesTmp_20160518/'
    head = 'T80Cam_T80_';  ext = '.tab'
    
    fname = []
    waveleng = []; transm = []
    D = []; N = []; Nb = []
    alpha = []; beta = []
    for i, filt in zip(np.arange(0,4), filset):
        fname.append(path + head + filt + ext)
        dummy, ddummy = np.loadtxt(fname[i], skiprows=1, unpack=True)
        waveleng.append(dummy)
        transm.append(ddummy) # Transmission curve of the 3 filters

        dummy = np.trapz(transm[i]*waveleng[i], dx=1)
        D.append(dummy)

        dummy = np.trapz(transm[i]*(waveleng[i]**2.), dx=1)
        N.append(dummy)

    N = np.array(N)
    D = np.array(D)
    alpha = N/D
    
    window = ((transm[0] > 0.) & (transm[1] > 0.)) #overlap window between nb and bbc filter
    lamb = np.median(waveleng[0][window])
    
    for i in np.arange(0,3):
        dummy = transm[i][window]*waveleng[i][window]
        Nb.append(np.median(dummy))

    Nb = np.array(Nb)
    beta = Nb/D

    return alpha, beta, lamb



#Computes the error on the continuum as estimated via
#the 3FM #(with both mocks3FM or jplus3FM). It assumes
#that alpha and beta coefficients have no error (could be super wrong)
def error3FM(ff, ee, aa, bb, ll):
    d = (aa[2]-aa[1])/(aa[0]-aa[2])
    dd = bb[0]*d + bb[1]
    sigma_l = (1./dd)*np.sqrt(ee[1]**2. + ee[2]**2. + (d**2.)*(ee[0]**2.+ee[2]**2.) )
    sigma_m = (1./(aa[0]-aa[1]) )*np.sqrt( ee[0]**2. + ee[1]**2. + (bb[0]*sigma_l)**2. )
    sigma_n = np.sqrt( (aa[2]*sigma_m)**2. + ee[2]**2. )
    sigma_c = np.sqrt( (ll*sigma_m)**2. + sigma_n**2. )

    return sigma_c*0.02




#Computes the continuum and line flux on jplus data given a 3-filters set
def jplus3FM(data, fset=setup['filters'], data_dir=setup['data']):
    
    alpha, beta, lamb = alphabeta(filset=fset)

    #mag to flux_lambda conversion
    flux_nb, dflux_nb = MtF(data[fset[0]][:,0], band=fset[0], dmags=data[fset[0]][:,1], unit='l')
    flux_bbc, dflux_bbc = MtF(data[fset[1]][:,0], band=fset[1], dmags=data[fset[1]][:,1], unit='l')
    flux_bbuc, dflux_bbuc = MtF(data[fset[2]][:,0], band=fset[2], dmags=data[fset[2]][:,1], unit='l')
    
    
    fline = ((flux_bbc-flux_bbuc) - ((alpha[1]-alpha[2])/(alpha[0]-alpha[2]))*(flux_nb-flux_bbuc))/((beta[0]*((alpha[2]-alpha[1])/(alpha[0]-alpha[2])) )+beta[1])
    M = (flux_nb - flux_bbuc - beta[0]*fline)/(alpha[0]-alpha[2])
    N = flux_bbuc - alpha[2]*M
    fcont = M*lamb + N

    fluxes = np.array([flux_nb, flux_bbc, flux_bbuc])
    errors = np.array([dflux_nb, dflux_bbc, dflux_bbuc])
    fcont_error = error3FM(fluxes, errors, alpha, beta, lamb)
    
    #in case some fline elements are nan, I need to mask
    flux_mask =  np.isnan(fline)
    fline = fline[~flux_mask]
    fcont = fcont[~flux_mask]
    flux_nb = flux_nb[~flux_mask]
    fcont_error = fcont_error[~flux_mask]

    return fline, fcont, flux_nb, fcont_error




#Computes the continuum and line flux on mocks data given a 3-filters set
def mocks3FM(data, fset=setup['filters'], data_dir=setup['data']):
        
    alpha, beta, lamb = alphabeta()
    
    #mag to flux_lambda conversion
    flux_nb = MtF(data[fset[0]][:,0], band=fset[0], unit='l')
    flux_bbc = MtF(data[fset[1]][:,0], band=fset[1], unit='l')
    flux_bbuc = MtF(data[fset[2]][:,0], band=fset[2], unit='l')
        
    fline = ((flux_bbc-flux_bbuc) - ((alpha[1]-alpha[2])/(alpha[0]-alpha[2]))*(flux_nb-flux_bbuc))/((beta[0]*((alpha[2]-alpha[1])/(alpha[0]-alpha[2])) )+beta[1])
    M = (flux_nb - flux_bbuc - beta[0]*fline)/(alpha[0]-alpha[2])
    N = flux_bbuc - alpha[2]*M
    fcont = M*lamb + N
    
    #in case some fline elements are nan, I need to mask
    flux_mask =  np.isnan(fline)
    fline = fline[~flux_mask]
    fcont = fcont[~flux_mask]
    flux_nb = flux_nb[~flux_mask]
    
    return fcont, fline, flux_nb


