
import os
import sys
import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt

import tools.settings
setup = tools.settings.set_up()


    

#Computes the continuum and line flux on jplus images given a 3-filters set
def threeFiltMeth(data):
    from threeFM import alphabeta
    
    alpha, beta, lamb = alphabeta(filset=setup['filters'])

    flux_nb = data[0]
    flux_bbc = data[1]
    flux_bbuc = data[2]
    
    
    fline = ((flux_bbc-flux_bbuc) - ((alpha[1]-alpha[2])/(alpha[0]-alpha[2]))*(flux_nb-flux_bbuc))/((beta[0]*((alpha[2]-alpha[1])/(alpha[0]-alpha[2])) )+beta[1])
    M = (flux_nb - flux_bbuc - beta[0]*fline)/(alpha[0]-alpha[2])
    N = flux_bbuc - alpha[2]*M
    fcont = M*lamb + N

    return fline, fcont


# This function reads 3 jplus images: 2 broad bands and 1 narrow band,
# performs a linear fit on the 2 broad bands and extract the NB excess.
# All this is done in a map-oriented way, so it constructs maps of NB excess from the images
# the input datafolder is just needed for saving the plot at the end
def read_jplus_image(inp_fit, zps, datafolder):
    from tools.jplus_filter_system import jplus_pivot
    
    cc = 2.999e18                               # in Angstrom
    mm = []

    plt.figure(num=None, figsize=( 11.692*2 , 8.267*1 )  )
    ax1 = plt.subplot2grid( ( 1 , 4 ) , (0, 0)  )
    ax2 = plt.subplot2grid( ( 1 , 4 ) , (0, 1), sharey=ax1 )
    ax3 = plt.subplot2grid( ( 1 , 4 ) , (0, 2), sharey=ax1 )
    ax4 = plt.subplot2grid( ( 1 , 4 ) , (0, 3), sharey=ax1 )

    for i in range(len(inp_fit)):
        if i == 0: accs = ax1 
        if i == 1: accs = ax2 
        if i == 2: accs = ax3 
        data, header = pf.getdata(inp_fit[i], 0, header=True)
        lpiv = jplus_pivot(band=setup['filters'][i])            # in Angstrom
        flu = data*(10**(-0.4*(zps[i] + 48.6)) )*(cc/(lpiv**2.))
        mm.append(flu)

        plot = accs.imshow(np.log10(flu), cmap='Blues')
        plt.colorbar(plot, ax=accs)
        accs.set_title(setup['filters'][i], fontsize=13)
        accs.set_xlabel('pixels', fontsize=13)
        accs.set_ylabel('pixels', fontsize=13)

    mm = np.array(mm)
    NBmap, contMap = threeFiltMeth(mm)
    
    plot = ax4.imshow(np.log10(NBmap), cmap='Reds')
    #ax4.imshow(np.log10(NBmap), cmap='Reds')
    #plt.colorbar(plot, ax4)
    ax4.set_title(setup['filters'][0]+' excess', fontsize=13)
    ax4.set_xlabel('pixels', fontsize=13)
    ax4.set_ylabel('pixels', fontsize=13)
    plt.tight_layout()
    plt.savefig(datafolder + 'filter_maps.png')
    #plt.show()
    plt.close()
    
    # SINGLE PLOT OF NB EXCESS
    plt.figure(figsize=(12,10))
    plt.imshow(np.log10(NBmap), cmap='Reds')
    plt.colorbar()
    plt.title(setup['filters'][0]+' excess')
    plt.savefig(datafolder + 'NBexcess_map.png')
    #plt.show()
    
    # print JPLUS header values
    # cp_nx = header['CRPIX1']   # number of central pixel in x dimension (ra dimension)
    # cp_ra = header['CRVAL1']   # ra of central pixel in x dim
    # cp_ny = header['CRPIX2']   # number of central pixel in y dim (dec dim)
    # cp_dec = header['CRVAL2']  # dec of central pixel in y dim
    # print cp_nx, cp_ra
    # print cp_ny, cp_dec

    sys.exit()
