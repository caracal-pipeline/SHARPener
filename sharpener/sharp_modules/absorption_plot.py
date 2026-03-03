__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys
import string
import os
import numpy as np
import yaml
import json
import glob
from astropy import wcs
from astropy.io import fits, ascii
from astropy import units as u
from astropy.table import Table, Column, MaskedColumn
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
# from mpdaf.obj import Spectrum, WaveCoord
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib.colors as mc


from astropy import units as u
from astropy.coordinates import Angle

import pypdf 
from sharpener.sharp_modules import convert_units as conv_units
from sharpener.sharp_modules import hi 

import logging

hi=hi.hi()

C = 2.99792458e5  # km/s
HI = 1.420405751e9  # Hz

####################################################################################################

def create_all_abs_plots(cfg_par):
    '''Function to create all absorption plots

    This function creates all absorption plots by calling abs_plot() for each source.
    Works as a kind of wrapper and is used by run_sharpener.py
    '''

    # get the list of spectra
    spectra = glob.glob(
        '{0:s}/*.txt'.format(cfg_par['general']['specdir']))
    # print(spectra)
    # check whether the previous step was successful
    if len(spectra) == 0:
        print("ERROR: No spectra found. Run spectrum extraction first")
        sys.exit(1)

    # go through all spectra
    #catalog_table = '{:s}{:s}'.format(cfg_par['general'].get('workdir'),cfg_par['source_catalog'].get('catalog_file'))
    # add sources
    mirCatalogFile = cfg_par['general']['absdir']+'mir_src_sharp.csv'
    catalog_table = '{:s}{:s}'.format(cfg_par['general'].get('absdir'),cfg_par['source_catalog'].get('catalog_file'))

    catalog_pybdsf = '{:s}{:s}'.format(cfg_par['general'].get('workdir'),cfg_par['source_catalog'].get('catalog_file'))
    
    #load the correct table accoridng to what has been used for catalog/sourcefinder
    if os.path.exists(catalog_pybdsf) and (cfg_par['source_catalog']['catalog']=='PYBDSF'):
        import Tigger
        from astropy.coordinates import Angle

        model = Tigger.load(catalog_pybdsf)
        sources = model.sources
        ra=[]
        dec=[]
        src_list=[]
        i=0
        for source in sources:
            ra_deg_angle  = Angle(np.rad2deg(source.pos.ra) * u.deg)
            dec_deg_angle = Angle(np.rad2deg(source.pos.dec) * u.deg)
            ra_hms = ra_deg_angle.to_string(unit=u.hourangle, sep=':')
            dec_dms = dec_deg_angle.to_string(unit=u.degree, sep=':')
            flux_cont.append(source.flux.I)

            src_list.append('{:s}_J{:s}{:s}{:s}.txt'.format(str(i),ra_hms.replace(':', ''),
                                                    '+' if source.pos.dec > 0.0 else '-',
                                                    dec_dms.replace(':', '')))
            i+=1

    elif os.path.exists(mirCatalogFile):
        sources = ascii.read(mirCatalogFile)
        src_list = []
        for i in range(len(sources)):
            src_list.append('{:s}_J{:s}.txt'.format(str(i),sources['J2000'][i]))

    elif os.path.exists(catalog_table) and (cfg_par['source_catalog']['catalog']=='NVSS'):
        sources = ascii.read(catalog_table)
        src_list =[]
        for i in range(len(sources)):
            src_list.append('{:s}_J{:s}.txt'.format(str(i),sources['NVSS'][i]))
    else:
        from astropy import units as u
        from astropy.coordinates import Angle
        catalog_table = '{:s}{:s}'.format(cfg_par['general'].get('workdir'),
                                                  cfg_par['source_catalog'].get('catalog_file'))
        vot = Table.read(catalog_table)
        i=0
        src_list=[]
        src_list_tmp=[]
        flux_cont=[]
        for row in vot:
            ra_deg_angle  = Angle((np.round(row['ra_peak'],2)) * u.deg)
            dec_deg_angle = Angle((np.round(row['dec_peak'],4)) * u.deg)
            ra_hms = ra_deg_angle.to_string(unit=u.hourangle, sep=':').split('.')[0]
            dec_dms = dec_deg_angle.to_string(unit=u.degree, sep=':').split('.')[0]

            J2000_name ='J{:s}{:s}'.format(ra_hms.replace(':', ''),dec_dms.replace(':', ''))
            src_list_tmp = '{:d}_{:s}.txt'.format(i,J2000_name)
            specName = cfg_par['general']['specdir']+src_list_tmp
            print(specName)
            if os.path.exists(specName):
                src_list.append(src_list_tmp)
                flux_cont.append(row['f_max'])

            i+=1
        # print(src_list)
        # src_list = os.path.basename(src_list)
    
    for i in range (0, len(src_list)):
        specName = cfg_par['general']['specdir']+os.path.basename(src_list[i])
        if os.path.exists(specName):
            abs_plot(specName, cfg_par)

    if cfg_par['abs_plot']['plot_contImage'] == True:
        plot_continuum(cfg_par)

    if cfg_par['abs_plot']['plot_detection_limits'] == True:
        plot_detection_limits(cfg_par,src_list,flux_cont)


def plot_detection_limits(cfg_par,src_list,flux_cont):
    '''Function to plot the detection limits given the continuum flux of the catalogue and the corresponding average noise in the spectra in sharpOut/spec
    '''
    nhi=np.zeros([len(src_list)],dtype=float)

    for i in range(0,len(src_list)):
        
        specName = cfg_par['general']['specdir']+os.path.basename(src_list[i])
        
        flux_cont_source = flux_cont[i]
       
        if os.path.exists(specName):
            spec_vec = ascii.read(specName)
            noise_spec =np.array(spec_vec[spec_vec.colnames[-1]][0],dtype=float)
            tau = hi.optical_depth(cfg_par['abs_plot']['sigma_detection_limit']*noise_spec,flux_cont_source)
            dv = cfg_par['abs_plot']['dv_detection_limit']
            nhi[i] = hi.nhi_abs(tau,dv)

    flux_cont = np.asarray(flux_cont)
    nhi = np.asarray(nhi)
    mask = (flux_cont > 0) & (nhi > 0) & (~np.isnan(flux_cont)) & (~np.isnan(nhi))
    x_data = np.log10(flux_cont[mask])
    y_data = np.log10(nhi[mask])

    # --- 2. Linear Regression (Log-Log Space) ---
    # p[0] is slope, p[1] is intercept
    p, cov = np.polyfit(x_data, y_data, 1, cov=True)
    slope, intercept = p
    
    # Calculate residuals to get the 1-sigma spread of the distribution
    fit_values = np.polyval(p, x_data)
    residual_std = np.std(y_data - fit_values)

    # Create a smooth line for the plot
    x_fit = np.log10(np.geomspace(flux_cont[mask].min(), flux_cont[mask].max(), 100))
    y_fit = np.polyval(p, x_fit)


    if cfg_par['abs_plot']['single_detections']:
        catalog_table = '{:s}{:s}'.format(cfg_par['general'].get('workdir'),
                                                  cfg_par['source_catalog'].get('catalog_file'))
        vot = Table.read(catalog_table)
        i=0
        flux_cont_det=[]
        for row in vot:
            flux_cont_det.append(row['f_max'])
        
        det_src = cfg_par['abs_plot']['single_detections']
        nhi_src=np.zeros([len(det_src)],dtype=float)
        flux_cont_det_src=np.zeros([len(det_src)],dtype=float)
        i=0
        for source in det_src:

            print(source)
            specName = cfg_par['general']['specdir']+source+'.txt'
            source_num = source.split('_')[0]
            flux_cont_det_src[i] = flux_cont_det[int(source_num)]

            if os.path.exists(specName):
                spec_vec = ascii.read(specName)

                print(specName)
                flux =np.array(spec_vec[spec_vec.colnames[1]],dtype=float)
                print(np.nanmin(flux),flux_cont_det[i])
                tau = hi.optical_depth(np.nanmin(flux),flux_cont_det_src[i])
                dv = 1.4
                nhi_src[i] = -hi.nhi_abs(tau,dv)
                print(nhi_src)
            i+=1
    flux_cont_det = np.asarray(flux_cont_det)
    nhi_det = np.asarray(nhi_src)
    print('#########')
    print(nhi_det)


    plt.rcParams.update({
            'font.family': 'serif',
            'text.usetex': True,
            'pgf.rcfonts': False,
            'text.latex.preamble': r'\usepackage{amsmath}\usepackage{amssymb}',
            'figure.facecolor': 'white',
            'xtick.direction': 'in',
            'ytick.direction': 'in',
            'xtick.top': True,
            'ytick.right': True,
            'axes.linewidth'      : 1.5,
            'lines.linewidth'     : 1.,
            'xtick.labelsize'     : 14,
            'ytick.labelsize'     : 14,
            'legend.fontsize'     : 10, 
            'xtick.direction'     :'in',
            'ytick.direction'     :'in',
            'xtick.major.size'    : 3,
            'xtick.major.width'   : 1.5,
            'xtick.minor.size'    : 2.5,
            'xtick.minor.width'   : 1.,
            'ytick.major.size'    : 3,
            'ytick.major.width'   : 1.5,
            'ytick.minor.size'    : 2.5,
            'ytick.minor.width'   : 1., 
        })

    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)

    # Scatter data (Lowered alpha to make the fit visible)
    ax.scatter(flux_cont, nhi, s=30, color='tab:red',marker='X', alpha=0.5, edgecolors='tab:red', label=r"$N(HI)^{abs}_{3\sigma,5km/s} = 1.9 \times 10^{18}\,T_s c_f\, \int \tau dv\, \mathrm{cm}^{-2}$")
    ax.scatter(flux_cont_det_src, nhi_det, marker='s', s=40, color='tab:blue', alpha=0.9, edgecolors='tab:blue', label=r"HI absorption detections [N(HI)$_{\mathrm peak, 1.4dv}$]")

    # Horizontal Line at 1.9e19
    ax.axhline(1.2e19, color='black', lw=1.5, ls='--', label=r'$N(HI)^{em}_{3\sigma,25km/s} = 1.9 \times 10^{19} \, \mathrm{cm}^{-2}$')

    # Plot the Erwin's 

    # ax.plot(10**x_fit, 10**y_fit, color='black', lw=2, label='Mean Linear Fit')

    # Plot the 1-Sigma Error Region
    # We add/subtract the residual_std in log space
    # ax.fill_between(10**x_fit, 10**(y_fit - residual_std), 10**(y_fit + residual_std), 
    #                 color='gray', alpha=0.3, label=r'$1\sigma$ distribution')

    # Styling
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$S_{\mathrm{c}} \text{ (Jy)}$', fontsize=16)
    ax.set_ylabel(r'$N_{\mathrm{HI}} \text{ (cm}^{-2}\text{)}$', fontsize=16)
    ax.legend(frameon=False, loc='upper right', fontsize=12)
    ax.minorticks_on()


    # 6. Saving logic
    outplot = "{0:s}{1:s}_abs_detection_limits.png".format(cfg_par['general']['plotdir'],cfg_par['general']['label'])
    if cfg_par['abs_plot']['plot_format'] == "pdf":
        plt.savefig(outplot.replace('.png', ".pdf"),
                    bbox_inches='tight')
    else:
        plt.savefig(outplot,
                    bbox_inches='tight', dpi=100)

    plt.close(fig)
def plot_continuum(cfg_par):
    '''Function to plot the continuum image from where spectra are extracted
    '''

    cont_im = cfg_par['general']['contname']
    if os.path.exists(cont_im):
        #load wcs system
        hdulist = fits.open(cont_im)  # read input
        # read data and header
        #what follows works for wcs, but can be written better
        prihdr = hdulist[0].header  
        w=wcs.WCS(prihdr)  

        # RS: the rest of the function requires only 2 axis images
        if w.naxis == 4:
            w = w.dropaxis(3)
            w = w.dropaxis(2)
            img = hdulist[0].data[0][0]
        elif w.naxis == 3:
            w = w.dropaxis(2)
   #         img = hdulist[0].data[0]    
            img = hdulist[0].data
    elif os.path.exists(cfg_par['general']['cubename']):
        
        cube_im = cfg_par['general']['cubename']
        #load wcs system
        hdulist = fits.open(cube_im)  # read input
        # read data and header
        #what follows works for wcs, but can be written better
        prihdr = hdulist[0].header  
        w=wcs.WCS(prihdr)  

        # RS: the rest of the function requires only 2 axis images
        if w.naxis == 4:
            w = w.dropaxis(3)
            w = w.dropaxis(2)
            img = hdulist[0].data[0][0]
            img = np.zeros(img.shape)

        elif w.naxis == 3:
            w = w.dropaxis(2)
            img = hdulist[0].data[0]
            img = np.zeros(img.shape)
    else:
        print("ERROR: No datacube found. Check configuration file")
        sys.exit(1)

    fig = plt.figure()
    ax = plt.subplot(projection=w)
    # ax.imshow(img, vmin=cfg_par[key]['clip'],
    #           vmax=np.max(img), norm=mc.LogNorm(cfg_par[key]['clip']), origin='lower')
    #ax.imshow(img, vmin=float(cfg_par[key]['clip'])/10000., vmax=np.min([np.max(img), float(cfg_par[key]['clip'])*1]), origin='lower')
    #fig = ax.imshow(img, vmin=0, vmax=float(cfg_par[key]['clip'])/5, origin = 'lower')

    figa = ax.imshow(img, norm=mc.SymLogNorm(float(cfg_par['source_finder']['clip'])/5.,
                                            vmin=float(cfg_par['source_finder']['clip'])/5.), origin='lower')

    # fig = ax.imshow(img, norm=mc.SymLogNorm(
    #   float(cfg_par[key]['clip'])*10), origin='lower')
    cbar = plt.colorbar(figa)
    cbar.set_label('Flux Density [Jy/beam]')
    #ax.imshow(img, vmin=, origin='lower')
    #ax.coords.grid(color='white', ls='solid')
    ax.coords[0].set_axislabel('Right Ascension')
    ax.coords[1].set_axislabel('Declination')
    ax.coords[0].set_major_formatter('hh:mm')
    ax.set_title("{0:s}".format(cfg_par['general']['workdir'].split('/')[-2]))
    #ax.coords[0].set_ticks(direction='in')
    #ax.coords[1].set_ticks(direction='in')
    # ax.tick_params(axis='both', bottom='on', top='on', left='on', right='on',
    #          which='major', direction='in')

    mirCatalogFile = cfg_par['general']['absdir']+'mir_src_sharp.csv'
    catalog_table = '{:s}{:s}'.format(cfg_par['general'].get('absdir'),cfg_par['source_catalog'].get('catalog_file'))
    catalog_pybdsf = '{:s}{:s}'.format(cfg_par['general'].get('workdir'),cfg_par['source_catalog'].get('catalog_file'))

    if os.path.exists(mirCatalogFile):
        src_list = ascii.read(mirCatalogFile)
        coord_list = SkyCoord(src_list['ra'], src_list['dec'], unit=(u.hourangle, u.deg), frame='fk5')

        for k in range(len(coord_list.ra)):
            ax.scatter(coord_list[k].ra.value, coord_list[k].dec.value, transform=ax.get_transform('fk5'),
                            edgecolor='red', facecolor='none')
            ax.annotate("{0:d}".format(k+1), xy=(coord_list[k].ra.value, coord_list[k].dec.value), xycoords=ax.get_transform('fk5'),
                                       xytext=(1, 1), textcoords='offset points', ha='left', color="white") 
  
    elif os.path.exists(catalog_pybdsf) and (cfg_par['source_catalog']['catalog']=='PYBDSF'):
        import Tigger
        from astropy.coordinates import Angle

        model = Tigger.load(catalog_pybdsf)
        sources = model.sources
        ra=[]
        dec=[]
        for source in sources:
            ra_deg_angle  = Angle(np.rad2deg(source.pos.ra) * u.deg)
            dec_deg_angle = Angle(np.rad2deg(source.pos.dec) * u.deg)
            ra.append(ra_deg_angle)
            dec.append(dec_deg_angle)

        coord_list = SkyCoord(ra,dec, unit=(u.deg, u.deg), frame='fk5')
 
        for k in range(len(coord_list.ra)):
            ax.scatter(coord_list[k].ra.value, coord_list[k].dec.value, transform=ax.get_transform('fk5'),
                            edgecolor='red', facecolor='none')
            ax.annotate("{0:d}".format(k+1), xy=(coord_list[k].ra.value, coord_list[k].dec.value), xycoords=ax.get_transform('fk5'),
                                       xytext=(1, 1), textcoords='offset points', ha='left', color="white")   

    elif os.path.exists(catalog_table) and (cfg_par['source_catalog']['catalog']=='NVSS'):

        src_list = ascii.read(catalog_table)
        ra = np.array(src_list['RAJ2000'],dtype=str)
        dec = np.array(src_list['DEJ2000'],dtype=str)
        pixels=np.zeros([len(ra),2])
        kk=[]
        
        cube_im = cfg_par['general']['cubename']
        #load wcs system
        hdulist = fits.open(cube_im)  # read input
        # read data and header
        #what follows works for wcs, but can be written better
        prihdr = hdulist[0].header  
        w=wcs.WCS(prihdr)  
        if w.naxis == 4:
            w = w.dropaxis(3)
            w = w.dropaxis(2)
        if w.naxis == 3:
            w = w.dropaxis(2)
        ra_vec=[]
        dec_vec=[]
        for i in xrange(0,len(ra)):
            if ra[i] == 'nan':
                pixels[i, 0]= np.nan
                pixels[i, 1]= np.nan
            else:
                ra_deg = conv_units.ra2deg(ra[i])
                dec_deg = conv_units.dec2deg(dec[i])
                px,py=w.wcs_world2pix(ra_deg,dec_deg,0)
                if (0 < round(px,0) < prihdr['NAXIS1'] and
                        0 < round(py,0) < prihdr['NAXIS2']):
                    kk.append(i)
                    ra_vec.append(ra[i])
                    dec_vec.append(dec[i])
                else:
                    pass

        coord_list = SkyCoord(ra_vec, dec_vec , unit=(u.hourangle, u.deg), frame='fk5')


        for k in range(len(coord_list.ra)):
            ax.scatter(coord_list[k].ra.value, coord_list[k].dec.value, transform=ax.get_transform('fk5'),
                            edgecolor='red', facecolor='none')
            ax.annotate("{0:s}".format(str(kk[k])), xy=(coord_list[k].ra.value, coord_list[k].dec.value), xycoords=ax.get_transform('fk5'),
                                       xytext=(1, 1), textcoords='offset points', ha='left', color="white") 

    output = "{0:s}{1:s}_continuum.png".format(cfg_par['general'].get(
        'plotdir'), cfg_par['general']['workdir'].split('/')[-2])

    if cfg_par['abs_plot']['plot_format'] == "pdf":
        fig.savefig(output.replace(".png", ".pdf"), bbox_inches='tight')
    else:
        fig.savefig(output, bbox_inches='tight', dpi=300)
    


def abs_plot(spec_name, cfg_par):
    '''
    Plots spectra of all radio sources found by find_src_imsad
    saved in basedir/beam/abs/spec.
    Plots are stored in basedir/beam/abs/plot

    IN
            Spectra extracted by spec_ex

    IN cfga
            abs_ex_plot_xaxis= ' '      #: X-axis units ['velocity','frequency']
            abs_ex_plot_yaxis= ' '      #: Y axis units ['flux','optical depth']
            #: plots line at redshift of source in spectrum redshift must be stored in table of load_src_csv
            abs_ex_plot_redsrc= True
            abs_ex_plot_title= True     #: plot title: J2000 name of radio source
            abs_ex_plot_format= ' '     #: format of plot ['.pdf','.jpeg','.png']

    OUT
            For each source outputs have the following name:
            J2000_xaxis-unit_yaxis-unit.plot_format = J220919.87+180920.17_vel_flux.pdf

    '''
    verb = cfg_par['general']['verbose']
    key = 'abs_plot'

    os.chdir(cfg_par['general']['specdir'])
    
    params = {
        'figure.autolayout' : True,
        'figure.facecolor': 'white',
        'pdf.fonttype'        : 3,
        # 'font.serif'          :'times',
        'font.style'          : 'normal',
        'font.weight'         : 'book',
        'font.size'           : 10,
        'axes.linewidth'      : 1.5,
        'lines.linewidth'     : 1.,
        'xtick.labelsize'     : 10,
        'ytick.labelsize'     : 10,
        'legend.fontsize'     : 10, 
        'xtick.direction'     :'in',
        'ytick.direction'     :'in',
        'xtick.major.size'    : 3,
        'xtick.major.width'   : 1.5,
        'xtick.minor.size'    : 2.5,
        'xtick.minor.width'   : 1.,
        'ytick.major.size'    : 3,
        'ytick.major.width'   : 1.5,
        'ytick.minor.size'    : 2.5,
        'ytick.minor.width'   : 1., 
        'text.usetex'         : True,
        'text.latex.preamble' : r'\usepackage{amsmath}',
        'text.latex.preamble' : r'\usepackage{lmodern}',    # latin modern, recommended to replace computer modern sans serif
        'text.latex.preamble' : r'\usepackage{helvet}',    # set the normal font here
         }
    plt.rcParams.update(params)

    #params = {
    #    'text.usetex': True,
    #}
   # rc('font', **{'family': 'serif', 'serif': ['serif']})

    # for i in xrange(0,len(np.atleast_1d(spec_src_name))):

    # load data and labels
    #	spec_name = spec_src_name[i]
    #	print spec_name
    if os.path.isfile(spec_name) == True:

        spec_vec = ascii.read(spec_name)
        x_data = np.array(spec_vec[spec_vec.colnames[0]], dtype=float)
        n_channels = np.size(x_data)

        # Set plot specs
        font_size = 16
        plt.ioff()
        plt.rc('xtick', labelsize=font_size-2)
        plt.rc('ytick', labelsize=font_size-2)

        fig, ax1 = plt.subplots(figsize=(12, 6))
        #fig = plt.figure(figsize=(9, 6))
        # fig.subplots_adjust(hspace=0.0)
        #gs = gridspec.GridSpec(1, 1)

        # Initialize subplots
        #ax1 = fig.add_subplot(gs[0])
        ax1.set_xlabel('')
        ax1.set_ylabel('')
        ax1.tick_params(axis='both', bottom='on', top='on',
                        left='on', right='on', which='major', direction='in')
        ax1.tick_params(axis='both', bottom='on', top='on',
                        left='on', right='on', which='minor', direction='in')

        flag_chans = cfg_par['spec_ex'].get('flag_chans', None)
        flag_chans1 = cfg_par['spec_ex'].get('flag_chans', None)

        if cfg_par['abs_plot'].get('zunit') == 'm/s':
            x_data /= 1e3
            ax1.set_xlabel(
                r'$cz\,(\mathrm{km}\,\mathrm{s}^{-1})$', fontsize=font_size)
        if cfg_par['abs_plot'].get('yunit') == 'tau':
            y_data = np.array(spec_vec[spec_vec.colnames[3]], dtype=float)
            y_sigma =np.array(spec_vec[spec_vec.colnames[4]])
            ylabh = ax1.set_ylabel(
            r'$\tau$', fontsize=font_size+2)
            ylabh.set_verticalalignment('center')
        else:
            y_data = np.array(spec_vec[spec_vec.colnames[1]], dtype=float)*1e3

            y_sigma = np.array(spec_vec[spec_vec.colnames[2]])*1e3
            ylabh = ax1.set_ylabel(
            r'$S_{\nu}$\,$[\mathrm{mJy}\,\mathrm{beam}^{-1}]$', fontsize=font_size)
            ylabh.set_verticalalignment('center')

        if cfg_par['abs_plot'].get('zunit') == 'MHz':
            x_data /= 1e6
            ax1.set_xlabel(r'Frequency [MHz]', fontsize=font_size)



        # Plot spectra
#                if self.abs_ex_plot_linestyle == 'step':
        # ax1.plot(x_data, y_data, color='black', linestyle='-')

        if flag_chans != None:
            flag_chans = np.array(flag_chans)
            if cfg_par['abs_plot'].get('zunit') == 'm/s':
                flag_chans = np.divide(flag_chans, 1e3)
            index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
            for k in xrange(1, len(flag_chans)):
                index_flags = (np.abs(x_data - flag_chans[k])).argmin()
                # y_data[index_flags] = 0.0
            y_data[index_flags_l:index_flags] = 0.0
        ax1.step(x_data, y_data, where='mid', color='black', linestyle='-',lw=1)

        # Calculate axis limits and aspect ratio
        x_min = np.nanmin(x_data)
        x_max = np.nanmax(x_data)
        y1_array = y_data[np.where((x_data > x_min) & (x_data < x_max))]
        if cfg_par[key]['fixed_scale']:
            y1_min = -5.*np.array(spec_vec[spec_vec.colnames[2]][0])*1e3*1.05
            y1_max = 5*np.array(spec_vec[spec_vec.colnames[2]][0])*1e3*1.05
        else:
            y1_min = np.nanmin(y_data)*1.1
            y1_max = np.nanmax(y_data)*1.1

        # Set axis limits
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(y1_min, y1_max)
        ax1.xaxis.labelpad = 6
        ax1.yaxis.labelpad = 10

        if flag_chans1 != None:
            flag_chans = np.array(flag_chans1)
            if cfg_par['abs_plot'].get('zunit') == 'm/s':
                flag_chans = np.divide(flag_chans, 1e3)
            index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
            for k in xrange(1, len(flag_chans)):
                index_flags = (np.abs(x_data - flag_chans[k])).argmin()
                ax1.fill_between([x_data[index_flags_l], x_data[index_flags]],
                                 y1_min, y1_max, facecolor='grey', alpha=0.3)

        # Plot noise
        ax1.fill_between(x_data, -y_sigma, y_sigma,
                         facecolor='grey', alpha=0.5,step='mid')

        # Plot stuff
        ax1.axhline(color='k', linestyle=':', zorder=0)
        if cfg_par[key]['fixed_scale']:
            ax1.axhline(color='k', linestyle=':', y=np.array(spec_vec[spec_vec.colnames[2]][0])*1e3)
            ax1.axhline(color='k', linestyle=':', y=-np.array(spec_vec[spec_vec.colnames[2]][0])*1e3)
            ax1.axhline(color='k', linestyle=':', y=-2*np.array(spec_vec[spec_vec.colnames[2]][0]*1e3))
            ax1.axhline(color='tab:green', linestyle=':', lw=1.5, y=-3*np.array(spec_vec[spec_vec.colnames[2]][0]*1e3))
            ax1.axhline(color='k', linestyle=':', y=-4*np.array(spec_vec[spec_vec.colnames[2]][0]*1e3))
            ax1.axhline(color='tab:red', linestyle=':', lw=1.5, y=-5*np.array(spec_vec[spec_vec.colnames[2]][0]*1e3))
            ax1.axvline(color='k', linestyle=':', x=x_data[int(len(x_data)/3)])
            ax1.axvline(color='k', linestyle=':', x=x_data[int(2*len(x_data)/3)])

        redshifts = cfg_par[key].get('redshift_sources', None)

        if len(redshifts) == 2:
            ax1.fill_between([redshifts[0], redshifts[1]], y1_min,
                             y1_max, facecolor='red', alpha=0.1)

        # Add title
        if cfg_par[key]['title'] == True:
            ax1.set_title(r"{0:s} (\#{1:d}): {2:s}".format(cfg_par['general']['label'], int(os.path.basename(spec_name).split('_')[0]), os.path.basename(spec_name).replace(
                '.txt', '').split('_')[-1]), fontsize=font_size+2)
            # if self.abs_ex_plot_title == True:
        #	ax1.set_title('%s' % (self.J2000_name[i]), fontsize=font_size+2)
        # ax1.axes.titlepad = 8

        # Add minor tick marks
        ax1.minorticks_on()

        # Save figure to file
        # name of plot is combination of beam number and name of extracted source
        outplot = os.path.basename(spec_name)
        # changed this to account for the different name convention
        # outplot = string.split(outplot,'.')[0]
        # outplot = cfg_par['general']['plotdir']+outplot+'.png'

        outplot = "{0:s}{1:s}_{2:s}".format(
            cfg_par['general']['plotdir'], cfg_par['general']['label'], outplot.replace('.txt', '_compact.png'))
        if cfg_par[key]['plot_format'] == "pdf":
            plt.savefig(outplot.replace('.png', ".pdf"),
                        bbox_inches='tight')
        else:
            plt.savefig(outplot,
                        bbox_inches='tight', dpi=100)

        plt.close("all")

        # also create multi-plot spectra
        if cfg_par[key]['detailed_plot']:
            # print(n_channels)

            # number of channels at which to split
            n_channel_per_plot = int(cfg_par[key]['channels_per_plot'])

            # get the number of plots
            n_plots = int(np.ceil(float(n_channels)/float(n_channel_per_plot)))

            # print(n_plots)

            n_rows = n_plots

            # add one row for the plot with full channel width
            fig, ax = plt.subplots(squeeze=False,
                ncols=1, nrows=n_rows, figsize=(10, 2*n_rows))
            fig.subplots_adjust(hspace=0.2)

            # ax1.annotate("Full spectrum", xy=(
            #     0.05, 0.95), xycoords='axes fraction', ha='left')

            # ax[1].annotate("Detailed spectrum", xy=(
            #     0.05, 0.95), xycoords='axes fraction', ha='left')

            if cfg_par[key]['title'] == True:
                ax[0][0].set_title(r"{0:s} (\#{1:d}): {2:s}".format(cfg_par['general']['label'], int(os.path.basename(spec_name).split('_')[0]), os.path.basename(spec_name).replace(
                    '.txt', '').split('_')[-1]), fontsize=font_size+2)

            # go through the rest of the plots and create them
            for plot_count in range(n_rows):

                # the chunk of data corresponding to the plot
                data_indices_min = plot_count * n_channel_per_plot
                data_indices_max = (plot_count+1) * n_channel_per_plot
                # print(data_indices_min)
                # print(data_indices_max)
                x_data_plot = x_data[data_indices_min:data_indices_max]
                y_data_plot = y_data[data_indices_min:data_indices_max]
                y_sigma_plot = y_sigma[data_indices_min:data_indices_max]

                # set the plot limits (only the x-axis needs to be adjusted)
                ax[plot_count][0].set_ylim(y1_min, y1_max)
                ax[plot_count][0].xaxis.labelpad = 6
                ax[plot_count][0].yaxis.labelpad = 10
                ax[plot_count][0].minorticks_on()
                ax[plot_count][0].tick_params(axis='both', bottom='on', top='on',
                                           left='on', right='on', which='major', direction='in')
                ax[plot_count][0].tick_params(axis='both', bottom='on', top='on',
                                           left='on', right='on', which='minor', direction='in')
                ylabh = ax[plot_count][0].set_ylabel(
                    r'S$_\nu$\,$[\mathrm{mJy}\,\mathrm{beam}^{-1}]$', fontsize=font_size)
                ylabh.set_verticalalignment('center')
                # adjust the plot range of the last plot to match the others if the number
                # of channels cannot be divided by the number of channels per plot without rest
                if plot_count == n_plots-1 and float(n_channels) % float(n_channel_per_plot) != 0:
                    # get channel spacing assuming the spacing is linear
                    channel_width = np.diff(x_data_plot)[0]

                    # number of missing channels
                    n_missing_channels = n_channel_per_plot - len(x_data_plot)

                    x_data_plot_min = np.min(
                        x_data_plot) + channel_width * n_missing_channels

                    # fill array up to the number of channels
                    # for k in range(n_missing_channels):
                    #     x_data_plot = np.append(
                    #         x_data_plot, x_data_plot[-1]+channel_width)
                    #     y_data_plot = np.append(y_data_plot, 0.)
                else:
                    x_data_plot_min = np.min(x_data_plot)

                x_data_plot_max = np.max(x_data_plot)
                #print(x_data_plot_max - x_data_plot_min)

                ax[plot_count][0].step(x_data_plot, y_data_plot, where='mid', color='black', linestyle='-')
                # Plot noise
                ax[plot_count][0].fill_between(x_data_plot, -y_sigma_plot, y_sigma_plot,
                         facecolor='grey', alpha=0.5,step='mid')
                ax[plot_count][0].set_xlim(x_data_plot_min, x_data_plot_max)

                # ax[plot_count].fill_between(x_data_plot, -y_sigma_plot, y_sigma_plot,
                #                             facecolor='grey', alpha=0.5)
                ax[plot_count][0].axhline(color='k', linestyle=':', zorder=0)

                # for the last plot add the x-axis label
                if plot_count == n_plots-1:

                    if cfg_par['spec_ex'].get('abs_plot') == 'm/s':
                        x_data /= 1e3
                        ax[plot_count][0].set_xlabel(
                            r'$cz\,(\mathrm{km}\,\mathrm{s}^{-1})$', fontsize=font_size)
                    elif cfg_par['spec_ex'].get('abs_plot') == 'MHz':
                        x_data /= 1e6
                        ax[plot_count][0].set_xlabel(
                            r'Frequency [MHz]', fontsize=font_size)

                    # in case the last plot contains less channels than the others
                    # if n_channels % n_channel_per_plot == 0:
                    #     n_channels_subplot = np.size(x_data_plot)
                    #     figwidth = fig.get_figwidth()
                    #     ax[plot_count].set_figwidth(
                    #         figwidth * n_channels_subplot/n_channel_per_plot)

            # name of plot
            outplot = os.path.basename(spec_name)
            outplot = "{0:s}{1:s}_{2:s}".format(
                cfg_par['general']['plotdir'], cfg_par['general']['label'], outplot.replace('.txt', '_detailed.png'))
            if cfg_par[key]['plot_format'] == "pdf":
                plt.savefig(outplot.replace('.png', ".pdf"),bbox_inches='tight')
            else:
                plt.savefig(outplot,bbox_inches='tight', dpi=100)

            plt.close("all")
        if verb == True:
            print('# Plotted spectrum of source ' + os.path.basename(spec_name)+'. #')
    else:
        print('# Missing spectrum of source ' + os.path.basename(spec_name)+'. #')


#######################################################################
##### Functions to plot spectra                                   #####
#######################################################################         

def plot_stack(cfg_par,stack_name):
    
    spec_vec=ascii.read(stack_name)

    x_data = np.array(spec_vec[spec_vec.colnames[0]], dtype=float)
    y_data = np.array(spec_vec[spec_vec.colnames[1]], dtype=float)
    y_sigma = np.array(spec_vec[spec_vec.colnames[2]], dtype=float)

    params = {
        'figure.autolayout' : True,
        'figure.facecolor': 'white',
        'pdf.fonttype'        : 3,
        # 'font.serif'          :'times',
        'font.style'          : 'normal',
        'font.weight'         : 'book',
        'font.size'           : 10,
        'axes.linewidth'      : 1.5,
        'lines.linewidth'     : 1.,
        'xtick.labelsize'     : 10,
        'ytick.labelsize'     : 10,
        'legend.fontsize'     : 10, 
        'xtick.direction'     :'in',
        'ytick.direction'     :'in',
        'xtick.major.size'    : 3,
        'xtick.major.width'   : 1.5,
        'xtick.minor.size'    : 2.5,
        'xtick.minor.width'   : 1.,
        'ytick.major.size'    : 3,
        'ytick.major.width'   : 1.5,
        'ytick.minor.size'    : 2.5,
        'ytick.minor.width'   : 1., 
        'text.usetex'         : True,
        'text.latex.preamble' : r'\usepackage{amsmath}',
        'text.latex.preamble' : r'\usepackage{lmodern}',    # latin modern, recommended to replace computer modern sans serif
        'text.latex.preamble' : r'\usepackage{helvet}',    # set the normal font here
         }
    plt.rcParams.update(params)

    
    line_size = 2

      # initialize figure
    font_size = 16
    plt.ioff()
    plt.rc('xtick', labelsize=font_size-2)
    plt.rc('ytick', labelsize=font_size-2)

    
    # Initialize subplots
    fig, ax1 = plt.subplots(figsize=(12, 6))

    ax1.set_xlabel(r'Velocity$\,[\mathrm{km}\,\mathrm{s}^{-1}]$', fontsize=params['font.size'])                
    
    # set y-label
    ylabh = ax1.set_ylabel(r'$\tau$', fontsize=params['font.size']+2)          
    ylabh.set_verticalalignment('center')

    # Calculate axis limits and aspect ratio
    x_min = np.min(x_data)
    x_max = np.max(x_data)
    y1_array = y_data[np.where((x_data>x_min) & (x_data<x_max))]
    y1_min = np.min(y1_array)*1.1
    y1_max = np.max(y1_array)*1.1

    # Set axis limits
    ax1.set_xlim(-cfg_par['stacking']['velrange']-20, cfg_par['stacking']['velrange']+20)
    ax1.set_ylim(y1_min, y1_max)
    ax1.xaxis.labelpad = 6
    ax1.yaxis.labelpad = 10

    # Plot spectra 
    # if abstack_plot_linestyle == 'step':
    ax1.step(x_data, y_data, where='mid', color='black', linestyle='-')
    # else:
    #     ax1.plot(x_data, y_data, color='black', linestyle='-')

    # Plot noise
    # ax1.fill_between(x_data, -y_sigma, y_sigma, facecolor='grey', alpha=0.5)

    #add vertical line at redshift of source
    ax1.axvline(color='k',linestyle=':', zorder = 0, lw=2)
    ax1.axhline(color='k', linestyle=':', zorder=0, lw=2)

    # Add title        
    if cfg_par['stacking']['plot_title'] != 'None':
        ax1.set_title(cfg_par['stacking']['plot_title'], fontsize=params['font.size']+2) 
    ax1.axes.titlepad = 8

    # Add minor tick marks
    ax1.minorticks_on()

    # Save figure to file
    out_stack_spec_plot= cfg_par['general']['plotdir']+'stacked_spectrum.'+cfg_par['abs_plot']['plot_format']   
    print(out_stack_spec_plot)
    plt.show()
    plt.savefig(out_stack_spec_plot,bbox_inches='tight', dpi=100)

    return 0
