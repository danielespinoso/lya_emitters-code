import sys
import matplotlib.pylab as plt
import matplotlib.colors as colors
import numpy as np

import tools

setup = tools.set_up()
filt = setup['filters']

violet = (0.6, 0.0, 0.6)
orange = (1.0, 0.5, 0.0)
yellow = (1.0, 1.0, 0.0)

def plot_colorMag_bytile(data_jpl, color, mask_selec=[], mask_tile=[],\
                         linea=[], mag_cut=20.5, tile_num=2525):
    
    sigline = linea(setup['marr'])
    plt.plot(data_jpl[filt[0]][mask_tile,0], color[mask_tile], 'og',\
             alpha=0.2, markersize=3, label='all data')   #plot all sources in tile 'tile_num'
    plt.plot(data_jpl[filt[0]][mask_selec,0], color[mask_selec], 'ob',\
             alpha=0.5, markersize=3, label='selected')   #plot selected sources in tile 'tile_num'
    plt.plot(setup['marr'], sigline, 'r--', linewidth=2)  #plot sigma_line
    plt.plot( (mag_cut, mag_cut), (-2., 15.), '--', c=violet, linewidth=2) #vertical line at mag_cut
    
    plt.axis([15.5, 24.5, -2, 7.])
    plt.text( 16., sigline[11]+0.25, r'$\Sigma$ = '+str(setup['sigma']),\
              color='r', fontweight='bold', fontsize=10)
    plt.text( (mag_cut+0.25), 6., 'S/N = '+str(setup['SN']),\
              color=violet, fontweight='bold', fontsize=10)
    plt.text( (mag_cut+0.25), 5.5, 'mag cut = '+str(int(mag_cut*100)/100.),\
              color=violet, fontweight='bold', fontsize=10)
    plt.title('jplus tile: '+str(tile_num))
    plt.legend(loc=2)
    plt.xlabel(str(filt[0])+'  [mags]')
    if setup['method'] == '3FM':
        plt.ylabel('['+str(filt[1])+'; '+str(filt[2])+'] - '+str(filt[0])+'   [mags]')
    elif setup['method'] == '2FM':
        plt.ylabel(str(filt[1])+' - '+str(filt[0])+'  [mags]')
    #plt.savefig(plotpath+'color-mag/'+setup['method']+'_'+filt[0]+'_color-mag_tile'+str(i)+'.png')
    plt.show()
    plt.close()


    

def plot_ColCol_bytile(data_jpl, data_mock=[], BroadBands=['gJAVA', 'rJAVA'],\
                       mask_sdss=[], mask_tile=[], mask_sel=[], tile_num=2525):

    bbc_color = data_jpl[BroadBands[0]][:,0] - data_jpl[filt[0]][:,0]
    bbuc_color = data_jpl[BroadBands[1]][:,0] - data_jpl[filt[0]][:,0]
    plt.plot(bbc_color[mask_tile], bbuc_color[mask_tile],'og', alpha=0.4, markersize=3, label='JPLUS')
    plt.plot(bbc_color[mask_sel], bbuc_color[mask_sel], 'ob', alpha=0.5,markersize=3,label='selected')

    if setup['plot_mock'] == True:
        mock_bbc_color = data_mock[BroadBands[0]][:,0] - data_mock[filt[0]][:,0]
        mock_bbuc_color = data_mock[BroadBands[1]][:,0] - data_mock[filt[0]][:,0]
        plt.plot(mock_bbc_color, mock_bbuc_color, 'or', alpha=0.4, markersize=3, label='mock')
    if setup['plot_sdssGal'] == True:
        GAL_bbc_color = data_jpl[BroadBands[0]][mask_sdss,0] - data_jpl[filt[0]][mask_sdss,0]
        GAL_bbuc_color = data_jpl[BroadBands[1]][mask_sdss,0] - data_jpl[filt[0]][mask_sdss,0]
        plt.plot(GAL_bbc_color, GAL_bbuc_color,'o',c=violet,alpha=0.3,markersize=3, label='sdss GAL')

    plt.axis([-0.5, 7., -0.5, 4.5])
    plt.plot( (setup['cCut'], setup['cCut']), (setup['cCut'], 4.5), '--', c=orange, linewidth=2)
    plt.plot( (setup['cCut'], 7.), (setup['cCut'], setup['cCut']), '--', c=orange, linewidth=2)
    plt.xlabel(filt[1]+' - '+filt[0]+'  [mags]')
    plt.ylabel(filt[2]+' - '+filt[0]+'  [mags]')
    plt.title('jplus tile: '+str(tile_num))
    plt.legend()
    #plt.savefig(plotpath+'color-color/'+setup['method']+'_'+filt[0]+'_UGcolor-color_tile'+str(i)+'.png')
    plt.show()
    plt.close()

    
    
def plot_CSTAR_bytile(data_jpl, mask_gaia=[], mask_sdss=[], mask_quasar=[],\
                       mask_tile=[], mask_sel=[], mag_cut=20.5, tile_num=2525):

    mask_sdss = ((mask_sdss) & (mask_tile))
    mask_quasar = ((mask_quasar) & (mask_tile))
    
    plt.axis([15, 30, -0.05, 1.1])
    #plt.plot(data_jpl[filt[0]][mask_tile,0], data_jpl['cstar'][mask_tile], 'og', markersize=4, alpha=0.1, label='jplus')
    plt.plot(data_jpl[filt[0]][mask_sel,0], data_jpl['cstar'][mask_sel], 'ob', markersize=4, alpha=0.6, label='selected')
    if setup['plot_gaia'] == True:
        plt.plot(data_jpl[filt[0]][mask_gaia,0], data_jpl['cstar'][mask_gaia], 'or', markersize=4, alpha=0.2, label='gaia')
    if setup['plot_sdssGal'] == True:
        plt.plot(data_jpl[filt[0]][mask_sdss,0], data_jpl['cstar'][mask_sdss], 'o', c=violet, markersize=4, alpha=0.3, label='sdss gal')
    if setup['plot_sdssQSO'] == True:
        plt.plot(data_jpl[filt[0]][mask_quasar,0], data_jpl['cstar'][mask_quasar], 'o', c=orange, markersize=4, alpha=0.4, label='sdss QSOs')
    #plt.plot( (mag_cut, mag_cut), (0, 10), 'c--', linewidth=2)
    #plt.text( (mag_cut+0.25), 9., 'mag cut = '+str(int(mag_cut*100)/100.),\
    #          color='c', fontweight='bold', fontsize=10)
    plt.legend()
    plt.xlabel(filt[0]+'  [mag]')
    plt.ylabel('CLASS_STAR')
    plt.title('jplus tile: '+str(tile_num))
    plt.show()
    plt.close()
    #sys.exit()

    
    
def plot_morpho_bytile(data_jpl, image_info, mask_gaia=[], mask_sdss=[], mask_quasar=[],\
                       mask_tile=[], mask_sel=[], mag_cut=20.5, tile_num=2525, PSF=False,\
                       MUMAX=True, pntLike_line=1.):
    # This function measures and plots the extended-ness of sources in a given tile
    # using two different methods: retriving it from the image PSF and from MU_MAX parameter
        
    mask_sdss = ((mask_sdss) & (mask_tile))
    mask_quasar = ((mask_quasar) & (mask_tile))
    mask_gaia = ((mask_gaia) & (mask_tile))
    
    mask_gq = ((mask_quasar) & (mask_gaia) & (mask_tile))
        
    fig = plt.figure(figsize=(10,8))

    # From PSF
    if PSF == True:
        psf_mask = ((image_info['tile_id'] == tile_num) & (image_info['filter_name'] == 'rJAVA'))
        morpho_pdf = data_jpl['fwhm'][:]/image_info['psf'][psf_mask][0]
        
        plt.axis([14, 22.5, 0, 10])
        
        plt.plot(data_jpl['rJAVA'][mask_tile,0], morpho_pdf[mask_tile], 'og',\
                 markersize=3, alpha=0.1, label='jplus')
        if setup['plot_gaia'] == True:
            plt.plot(data_jpl['rJAVA'][mask_gaia,0], morpho_pdf[mask_gaia], 'or',\
                     markersize=4, alpha=0.1, label='gaia')
        if setup['plot_sdssGal'] == True:
            plt.plot(data_jpl['rJAVA'][mask_sdss,0], morpho_pdf[mask_sdss], 'o',\
                     c=violet, markersize=4, alpha=0.35, label='sdss gal')
        if setup['plot_sdssQSO'] == True:
            plt.plot(data_jpl['rJAVA'][mask_quasar,0], morpho_pdf[mask_quasar], 'o',\
                     color=yellow, markersize=4, alpha=0.6, label='sdss QSOs')
        plt.plot(data_jpl['rJAVA'][mask_sel,0], morpho_pdf[mask_sel], 'ob', markersize=4,\
                 alpha=0.8, label='selected')
        
        plt.plot( (mag_cut, mag_cut), (0, 10), 'k--', linewidth=2)
        plt.text( (mag_cut-1.6), 9., 'mag cut = '+str(int(mag_cut*100)/100.),\
                  color='k', fontweight='bold', fontsize=10)
        plt.legend(loc=2)
        plt.xlabel('rJAVA'+'  [mag]')
        plt.ylabel('fwhm / PSF')
        plt.title('jplus tile: '+str(tile_num))
        #plt.savefig(setup['plots']+'FWHM_PSF_vs_'+filt[0]+'_tile'+str(tile_num)+'.png')
        plt.show()
        plt.close()
    
    elif MUMAX == True:
        # From MU_MAX parameter (mag/arcsec^2 of the bightest pixel)
        pxl = 0.55*0.55  # jplus pixel area in arcsec^2
        mumax = data_jpl['mu_max_r'] - 2.5*np.log10(pxl)
        extdness = mumax - data_jpl['rJAVA'][:,0]

        on_x = extdness[mask_tile]
        on_y = data_jpl['rJAVA'][mask_tile,0]

        H , y, x = np.histogram2d(on_x, on_y, bins = 150, normed = True)
        plt.pcolor(x, y, H, norm=colors.LogNorm(), cmap='Greens')
        
        if setup['plot_gaia'] == True:
            plt.plot(data_jpl['rJAVA'][mask_gaia,0], extdness[mask_gaia], 'or',\
                     markersize=2, alpha=0.1, label='gaia')
        if setup['plot_sdssGal'] == True:
            plt.plot(data_jpl['rJAVA'][mask_sdss,0], extdness[mask_sdss], 'o',\
                     c=violet, markersize=4, markeredgecolor='k', alpha=0.4, label='sdss gal')
        if setup['plot_sdssQSO'] == True:
            plt.plot(data_jpl['rJAVA'][mask_quasar,0], extdness[mask_quasar], 'o',\
                     color=yellow, markersize=4, markeredgecolor='k', alpha=0.7, label='sdss QSOs')

        # PLOT JPLUS OBJECTS IN BOTH GAIA AND SDSS_QSO
        # plt.plot(data_jpl['rJAVA'][mask_gq,0], extdness[mask_gq], 'o', color=orange,\
        #         markersize=2, alpha=0.6, label='gaia & QSOs')
            
        plt.plot(data_jpl['rJAVA'][mask_sel,0], extdness[mask_sel], 'ob',\
                 markersize=4, alpha=0.8, label='selected')

        plt.plot((12., 26.), (pntLike_line, pntLike_line), 'r--', linewidth=2)
        plt.plot((12., 26.), (pntLike_line+0.2, pntLike_line+0.2), 'k--', linewidth=1)
        # plt.plot((12., 26.), (pntLike_line-0.2, pntLike_line-0.2), 'k--', linewidth=1)
        
        plt.colorbar(label = '"Extended-ness"')
        plt.title('jplus tile: '+str(tile_num))
        plt.xlabel('rJAVA'+'  [mag]')
        plt.ylabel('MU_MAX(rJAVA)  -  rJAVA'+'  [mag]')
        plt.legend(loc=1)
        plt.savefig(setup['plots']+'MUMAX_tile'+str(tile_num)+'compactness_line.png')
        #plt.show()
        plt.close()
