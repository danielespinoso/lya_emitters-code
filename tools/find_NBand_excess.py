import numpy as np
import sys

from tools.converters import magtoflux as MtF

import tools.settings
setup = tools.settings.set_up()

 
def find_lineExcess(data, fset, sigma=setup['line_ex_sigma']):
    cln = ( (data[fset[0]][:,0] > 0) & (data[fset[1]][:,0] > 0) )
    fl_nb,dfl_nb=MtF(data[fset[0]][:,0], band=fset[0], dmags=data[fset[0]][:,1], unit='l')
    fl_bb_lft,dfl_lft=MtF(data[fset[1]][:,0],band=fset[1],dmags=data[fset[1]][:,1],unit='l')
    fl_bb_rgt,dfl_rgt=MtF(data[fset[2]][:,0],band=fset[2],dmags=data[fset[2]][:,1],unit='l')
    
    cond1 = (fl_nb > fl_bb_lft + sigma*dfl_lft)
    cond2 = (fl_nb > fl_bb_lft + sigma*dfl_nb)
    cond3 = (fl_nb > fl_bb_rgt + sigma*dfl_rgt)
    cond4 = (fl_nb > fl_bb_rgt + sigma*dfl_nb)
    true_excess = ( (cond1) & (cond2) & (cond3) & (cond4) )


    mask = ((true_excess) & (cln))
    
    # import matplotlib.pyplot as plt
    # for i in np.arange(0,len(data[fset[0]][:,0])):
    #     if mask[i] == True:
    #         for j in np.arange(0, 12):
    #             baan = tools.jpflt(j)
    #             lamb = tools.jplus_pivot(band=baan)
    #             flu, dflu = tools.singleMtoF(data[baan][i,0],band=baan,dmags=data[baan][i,1],unit='l')
    #             plt.errorbar(lamb, flu, yerr=dflu, c='k', fmt='s', markersize=4)
                
    #         baan = fset[0]
    #         lamb = tools.jplus_pivot(band=baan)
    #         flu, dflu = tools.singleMtoF(data[baan][i,0],band=baan,dmags=data[baan][i,1],unit='l')
    #         plt.errorbar(lamb, flu, yerr=dflu, c='r', fmt='s', markersize=4)
    #         plt.title(' "'+str(mask[i])+'" excess in NB: '+fset[0])
    #         plt.show()
    #         plt.close()

    return mask
