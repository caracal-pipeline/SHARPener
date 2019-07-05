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
from mpdaf.obj import Spectrum, WaveCoord
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib.colors as mc

import logging


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
        '{0:s}*.txt'.format(cfg_par['general']['specdir']))

    # check whether the previous step was successful
    if len(spectra) == 0:
        print("ERROR: No spectra found. Abort")
        sys.exit(1)

    # go through all spectra
    catalog_table = '{:s}{:s}'.format(cfg_par['general'].get('workdir'),cfg_par['source_catalog'].get('catalog_file'))
    if os.path.exists(catalog_table) and os.path.exists(cfg_par['general']['contname']):
        plot_continuum(cfg_par)

    for i in range(len(spectra)):
        spectra[i] = os.path.basename(spectra[i])
        abs_plot(spectra[i], cfg_par)


def plot_continuum(cfg_par):
    '''Function to plot the continuum image from where spectra are extracted
    '''

    cont_im = os.path.basename(cfg_par['general']['contname'])
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
        img = hdulist[0].data[0]



    ax = plt.subplot(projection=w)
    # ax.imshow(img, vmin=cfg_par[key]['clip'],
    #           vmax=np.max(img), norm=mc.LogNorm(cfg_par[key]['clip']), origin='lower')
    #ax.imshow(img, vmin=float(cfg_par[key]['clip'])/10000., vmax=np.min([np.max(img), float(cfg_par[key]['clip'])*1]), origin='lower')
    #fig = ax.imshow(img, vmin=0, vmax=float(cfg_par[key]['clip'])/5, origin = 'lower')

    fig = ax.imshow(img, norm=mc.SymLogNorm(float(cfg_par['source_finder']['clip'])/5.,
                                            vmin=float(cfg_par['source_finder']['clip'])/5.), origin='lower')

    # fig = ax.imshow(img, norm=mc.SymLogNorm(
    #   float(cfg_par[key]['clip'])*10), origin='lower')
    cbar = plt.colorbar(fig)
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

    # add sources
    mirCatalogFile = cfg_par['general']['absdir']+'mir_src_sharpener.csv'
    catalog_table = '{:s}{:s}'.format(cfg_par['general'].get('workdir'),cfg_par['source_catalog'].get('catalog_file'))

    if os.path.exists(mirCatalogFile):
        src_list = ascii.read(mirCatalogFile)
        coord_list = SkyCoord(src_list['ra'], src_list['dec'], unit=(u.hourangle, u.deg), frame='fk5')

    
    elif os.path.exists(catalog_table) and (cfg_par['source_catalog']['enable']==True):
        import Tigger
        from astropy.coordinates import Angle

        model = Tigger.load(catalog_table)
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

    output = "{0:s}{1:s}_continuum.png".format(cfg_par['general'].get(
        'plotdir'), cfg_par['general']['workdir'].split('/')[-2])

    if cfg_par['general']['plot_format'] == "pdf":
        plt.savefig(output.replace(".png", ".pdf"), overwrite=True, bbox_inches='tight')
    else:
        plt.savefig(output, overwrite=True, bbox_inches='tight', dpi=300)
    


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
        'text.usetex': True,
        'text.latex.unicode': True
    }
    rc('font', **{'family': 'serif', 'serif': ['serif']})
    plt.rcParams.update(params)

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

        fig, ax1 = plt.subplots(figsize=(10, 5))
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

        if cfg_par['spec_ex'].get('zunit') == 'm/s':
            x_data /= 1e3
            ax1.set_xlabel(
                r'$cz\,(\mathrm{km}\,\mathrm{s}^{-1})$', fontsize=font_size)
        y_data = np.array(spec_vec[spec_vec.colnames[1]], dtype=float)*1e3
        y_sigma = np.array(spec_vec[spec_vec.colnames[2]])*1e3

        if cfg_par['spec_ex'].get('zunit') == 'MHz':
            x_data /= 1e6
            ax1.set_xlabel(r'Frequency [MHz]', fontsize=font_size)

        ylabh = ax1.set_ylabel(
            r'S\,$[\mathrm{mJy}\,\mathrm{beam}^{-1}]$', fontsize=font_size)
        ylabh.set_verticalalignment('center')

        # Plot spectra
#                if self.abs_ex_plot_linestyle == 'step':
        # ax1.plot(x_data, y_data, color='black', linestyle='-')

        if flag_chans != None:
            flag_chans = np.array(flag_chans)
            if cfg_par['spec_ex'].get('zunit') == 'm/s':
                flag_chans = np.divide(flag_chans, 1e3)
            index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
            for k in xrange(1, len(flag_chans)):
                index_flags = (np.abs(x_data - flag_chans[k])).argmin()
                # y_data[index_flags] = 0.0
            y_data[index_flags_l:index_flags] = 0.0
        ax1.step(x_data, y_data, where='mid', color='black', linestyle='-')

        # Calculate axis limits and aspect ratio
        x_min = np.min(x_data)
        x_max = np.max(x_data)
        y1_array = y_data[np.where((x_data > x_min) & (x_data < x_max))]
        if cfg_par[key]['fixed_scale']:
            y1_min = -50
            y1_max = 50
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
            if cfg_par['spec_ex'].get('zunit') == 'm/s':
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

        redshifts = cfg_par[key].get('redshift_sources', None)
        if len(redshifts) == 2:
            ax1.fill_between([redshifts[0], redshifts[1]], y1_min,
                             y1_max, facecolor='red', alpha=0.1)

        # Add title
        if cfg_par[key]['title'] == True:
            ax1.set_title("{0:s} (\#{1:d}): {2:s}".format(cfg_par['general']['label'], int(os.path.basename(spec_name).split('_')[0]), os.path.basename(spec_name).replace(
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
        if cfg_par['general']['plot_format'] == "pdf":
            plt.savefig(outplot.replace('.png', ".pdf"),
                        overwrite=True, bbox_inches='tight')
        else:
            plt.savefig(outplot,
                        overwrite=True, bbox_inches='tight', dpi=100)

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
                ax[0][0].set_title("{0:s} (\#{1:d}): {2:s}".format(cfg_par['general']['label'], int(os.path.basename(spec_name).split('_')[0]), os.path.basename(spec_name).replace(
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
                    r'S\,$[\mathrm{mJy}\,\mathrm{beam}^{-1}]$', fontsize=font_size)
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

                    if cfg_par['spec_ex'].get('zunit') == 'm/s':
                        x_data /= 1e3
                        ax[plot_count][0].set_xlabel(
                            r'$cz\,(\mathrm{km}\,\mathrm{s}^{-1})$', fontsize=font_size)
                    elif cfg_par['spec_ex'].get('zunit') == 'MHz':
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
            if cfg_par['general']['plot_format'] == "pdf":
                plt.savefig(outplot.replace('.png', ".pdf"),
                            overwrite=True, bbox_inches='tight')
            else:
                plt.savefig(outplot,
                            overwrite=True, bbox_inches='tight', dpi=100)

            plt.close("all")
        if verb == True:
            print('# Plotted spectrum of source ' + os.path.basename(spec_name)+'. #')
    else:
        print('# Missing spectrum of source ' + os.path.basename(spec_name)+'. #')
