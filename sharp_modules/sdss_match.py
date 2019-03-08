#!/usr/bin/python2.7

from astropy.nddata.utils import Cutout2D
import matplotlib.colors as mc
import matplotlib.pyplot as plt
""" Module to get SDSS sources 

This module will take an Apertif image and determine the SDSS
sources that fit into the coordinate range and the frequency
range (i.e., redshift range). The remaining SDSS catalogue
entries will be saved or simply printed. If requested an image
will be plotted with the SDSS and radio sources in the field.

Issue: does not work in parrallel mode

__author__ = "Robert Schulz"
__copyright__ = "Robert Schulz"
__email__ = "schulz@astron.nl"

"""

# import modules
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack, hstack, Column
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.sdss import SDSS
import argparse
import sys
import os
import time
import matplotlib
matplotlib.use('Agg')


def get_sdss_sources(cfg_par):
    """This function extracts SDSS sources

    Notes:
        The function extracts SDSS sources from the SDSS catalogue
        based on their location in the Apertif image. It uses the cube
        as a reference rather than the continuum as the former is the
        reference for extracting the spectrum. Sometimes the continuum image
        is larger as required by cleaning.
        The resulting SDSS catalogue will be saved and can be matched to the
        radio sources depending on the settings for this module.
        It can also create a plot with all sources in it. The plot is cropped
        to the size of a channel image from the cube assuming they have both the
        same center

    Parameters:
        cfg_par : settings from the config file
    """

    print("# Running SDSS source finder")

    key = 'sdss_match'

    # HI frequency
    freq_hi = 1.420405751e9 * u.Hz

    # Getting directories and convert files if necessary
    # ++++++++++++++++++++++++++++++++++++++++++++++++++
    workdir = cfg_par['general']['workdir']
    os.chdir(workdir)
    cont_im = os.path.basename(cfg_par['general']['contname'])
    cube_im = os.path.basename(cfg_par['general']['cubename'])

    # cannot use cfg_par, probably because file name would be too long for miriad
    output_file_name = "{0:s}abs/{1:s}_sdss_src.csv".format(
        workdir, workdir.split('/')[-2])

    # Setting for redshift range from parameter file
    # ++++++++++++++++++++++++++++++++++++++++++++++
    freq_min = cfg_par[key]['freq_min'] * u.megaHertz
    freq_max = cfg_par[key]['freq_max'] * u.megaHertz

    # get corresponding redshift range
    z_max = freq_hi / freq_min - 1.
    z_min = freq_hi / freq_max - 1.
    if z_min < 0:
        z_min = 0

    print(
        "Searching in frequency range: {0:.3g} - {1:.3g}".format(freq_min, freq_max))
    print(
        "Corresponding to redshift range: {0:.3f} - {1:.3f}".format(z_min, z_max))

    # Process cube
    # ++++++++++++
    # better than continuum as it has all the relevant information and the right size

    # open fits file
    fits_hdulist = fits.open(cube_im)

    # get WCS header
    wcs = WCS(fits_hdulist[0].header)

    # if the number of axis are less than 4, the script will not work
    if wcs.naxis < 4:
        print("ERROR: The script requires a fits image with 4 axis. Abort")
        return -1

    # size of image
    image_shape = np.shape(fits_hdulist[0].data)

    # get image
    image = fits_hdulist[0].data[0][0]

    # get coordinates from image
    image_coordinates_1 = wcs.wcs_pix2world([[0, 0, 0, 0]], 1)
    image_coordinates_2 = wcs.wcs_pix2world(
        [[int(image_shape[-1]/2), int(image_shape[-2]/2), 0, 0]], 1)
    image_coordinates_3 = wcs.wcs_pix2world(
        [[image_shape[-1]-1, image_shape[-2]-1, 0, 0]], 1)

    # convert minimum and maximum RA and DEC values to SkyCoordinates
    image_coordinates_min = SkyCoord(
        image_coordinates_1[0][0], image_coordinates_1[0][1], unit='deg')
    image_coordinates_center = SkyCoord(
        image_coordinates_2[0][0], image_coordinates_2[0][1], unit='deg')
    image_coordinates_max = SkyCoord(
        image_coordinates_3[0][0], image_coordinates_3[0][1], unit='deg')

    # get separation
    image_radius = image_coordinates_max.separation(image_coordinates_min)
    print("Maximum angular radius of image: {0:.5f}deg".format(
        image_radius.value))

    # Get SDSS sources for image
    # ++++++++++++++++++++++++++
    print("Query SDSS catalogue")

    start_time_query = time.time()
    # because a radius is used, this step is not enough
    # and the SDSS have to further processed
    sdss_cat = SDSS.query_region(
        image_coordinates_center, radius=image_radius, spectro=True, timeout=int(cfg_par[key]['sdss_query_timeout']))

    print("Query SDSS catalogue ... Done {0:.1f}s".format(
        time.time()-start_time_query))

    # check that all sources are within the coordinate range of the image
    print("Matching SDSS sources to image size")
    sdss_cat = sdss_cat[np.where((sdss_cat['ra'] > image_coordinates_max.ra.value) & (
        sdss_cat['ra'] < image_coordinates_min.ra.value))]
    sdss_cat = sdss_cat[np.where((sdss_cat['dec'] > image_coordinates_min.dec.value) & (
        sdss_cat['dec'] < image_coordinates_max.dec.value))]

    # get the number of sources
    n_sdss_src = np.size(sdss_cat['ra'])

    if n_sdss_src == 0:
        print("# NO SDSS sources found. Abort")
        return -1

    # Process SDSS sources
    # ++++++++++++++++++++++++++
    print("## Process SDSS catalogue")

    # Match the SDSS sources redshift to the redshift range
    print("Matching SDSS sources to redshift range")
    sdss_cat = sdss_cat[np.where(
        (sdss_cat['z'] > z_min) & (sdss_cat['z'] < z_max))]

    # check that all sources are within the coordinate range of the image
    print("Matching SDSS sources to image size")
    sdss_cat = sdss_cat[np.where((sdss_cat['ra'] > image_coordinates_max.ra.value) & (
        sdss_cat['ra'] < image_coordinates_min.ra.value))]
    sdss_cat = sdss_cat[np.where((sdss_cat['dec'] > image_coordinates_min.dec.value) & (
        sdss_cat['dec'] < image_coordinates_max.dec.value))]

    # get the number of sources
    n_sdss_src = np.size(sdss_cat['ra'])

    # Removing those sources where there are no image values
    # Necessary, because of the shape of the mosaic
    print("Removing SDSS sources from where images is empty")
    src_list_keep = np.array([], dtype=int)
    for k in range(n_sdss_src):

        # get pixel from coordinates
        image_pixel = wcs.wcs_world2pix(
            [[sdss_cat['ra'][k], sdss_cat['dec'][k], 0, 0]], 1)

        if image_pixel[0][1] < image_shape[-2] and image_pixel[0][0] < image_shape[-1]:

            # get the image value
            image_value = image[int(
                image_pixel[0][1])][int(image_pixel[0][0])]

            if not np.isnan(image_value):
                src_list_keep = np.append(src_list_keep, k)

    sdss_cat = sdss_cat[src_list_keep]

    # get the number of sources
    n_sdss_src = np.size(sdss_cat['ra'])

    if n_sdss_src == 0:
        print('No SDSS source found.')
        sys.exit(1)
    else:
        print("There are {0:d} SDSS source in the given image field and redshift range".format(
            n_sdss_src))

    # save file or just print it
    sdss_cat.write(output_file_name, format="csv", overwrite=True)

    # close file
    fits_hdulist.close()

    # match sdss and radio
    if cfg_par[key]['match_cat']:

        print("Matching radio and SDSS catalogues")
        match_sdss_to_radio(cfg_par)

    # plot sdss sources
    if cfg_par[key]['plot_image']:

        print("Plotting SDSS and radio")

        # open continuum image
        # ++++++++++++++++++++
        fits_hdulist = fits.open(cont_im)

        # get WCS header
        wcs = WCS(fits_hdulist[0].header)

        # rest of the function requires only 2 axis images
        if wcs.naxis == 4:
            wcs = wcs.dropaxis(3)
            wcs = wcs.dropaxis(2)
            img = fits_hdulist[0].data[0][0]
        elif wcs.naxis == 3:
            wcs = wcs.dropaxis(2)
            img = fits_hdulist[0].data[0]
        else:
            img = fits_hdulist[0].data

        img_shape = img.shape

        # Crop image using cube
        # +++++++++++++++++++++
        fits_hdulist_cube = fits.open(cube_im)

        # get WCS header of cube
        wcs_cube = WCS(fits_hdulist_cube[0].header)

        # get shape of first image in cube
        if wcs_cube.naxis == 4:
            img_ch1_shape = fits_hdulist_cube[0].data[0][0].shape
        else:
            img_ch1_shape = fits_hdulist_cube[0].data[0].shape

        # cutout only if cube is smaller
        if img_ch1_shape[0] < img_shape[0] and img_ch1_shape[1] < img_shape[1]:

            # center of continuum image which is assumed to be the same as center of cube
            position = SkyCoord(
                wcs_cube.wcs.crval[0] * u.deg, wcs_cube.wcs.crval[1] * u.deg, frame='fk5')

            # create the cutout
            cutout = Cutout2D(img, position, img_ch1_shape, wcs=wcs)

            # overwrite existing wcs and image with new one from cutout
            wcs = cutout.wcs
            img = cutout.data

        # set up plot
        ax = plt.subplot(projection=wcs)

        #fig = ax.imshow(img, vmin=0, vmax=float(cfg_par['source_finder']['clip'])/5, origin = 'lower')

        fig = ax.imshow(img, norm=mc.SymLogNorm(float(cfg_par['source_finder']['clip'])/5.,
                                                vmin=float(cfg_par['source_finder']['clip'])/5.), origin='lower')

        cbar = plt.colorbar(fig)
        cbar.set_label('Flux Density [Jy/beam]')

        ax.coords[0].set_axislabel('Right Ascension')
        ax.coords[1].set_axislabel('Declination')
        ax.coords[0].set_major_formatter('hh:mm')
        ax.set_title("{0:s}".format(
            cfg_par['general']['workdir'].split('/')[-2]))

        # read the radio sources
        radio_src = Table.read("abs/mir_src_sharpener.csv")

        # number of radio sources
        n_radio_src = np.size(radio_src['ra'])

        # go through the sources to plot them
        print("Plotting radio sources")
        # create a SkyCoord object to convert the coordinates from string to float
        coord_radio_src = SkyCoord(
            radio_src['ra'], radio_src['dec'], unit=(u.hourangle, u.deg), frame='fk5')

        # go through radio sources and plot them
        for k in range(n_radio_src):

            # check that the radio sources are in the plot range
            if (coord_radio_src[k].ra.value >= image_coordinates_max.ra.value) and (coord_radio_src[k].ra.value <= image_coordinates_min.ra.value) and (coord_radio_src[k].dec.value >= image_coordinates_min.dec.value) and (coord_radio_src[k].dec.value <= image_coordinates_max.dec.value):

                # plot the sources
                ax.plot(coord_radio_src[k].ra.value, coord_radio_src[k].dec.value, transform=ax.get_transform('fk5'),
                        marker='o', markeredgecolor='red', markerfacecolor='none', zorder=3, markersize=5)

                # add the source number
                ax.annotate("{0:d}".format(k+1), xy=(coord_radio_src[k].ra.value, coord_radio_src[k].dec.value), xycoords=ax.get_transform('fk5'),
                            xytext=(1, 1), textcoords='offset points', ha='left', color="white")

        print("Plotting sdss sources")
        # add the SDSS sources

        # create a SkyCoord object to convert the coordinates from string to float
        coord_sdss_src = SkyCoord(
            sdss_cat['ra'], sdss_cat['dec'], unit=(u.deg, u.deg), frame='fk5')

        # plot the sdss sources
        for k in range(n_sdss_src):

            ax.plot(coord_sdss_src[k].ra.value, coord_sdss_src[k].dec.value, transform=ax.get_transform('fk5'), marker="s",
                    markeredgecolor='lightgray', markerfacecolor='none', markersize=4, zorder=1)

        # check if there was a match performed with SDSS and radio
        if cfg_par[key]['match_cat']:
            # file name
            radio_sdss_src_cat_file = "{0:s}abs/{1:s}_radio_sdss_src.csv".format(
                workdir, workdir.split('/')[-2])

            radio_sdss_src_cat = Table.read(
                radio_sdss_src_cat_file, format="csv")

            # get indices for which there was a match
            match_indices = np.where(radio_sdss_src_cat['sdss_ra'] != 0.)[0]

            for index in match_indices:
                # create a SkyCoord object to convert the coordinates from string to float
                coord_sdss_src = SkyCoord(
                    radio_sdss_src_cat[index]['sdss_ra'], radio_sdss_src_cat[index]['sdss_dec'], unit=(u.deg, u.deg), frame='fk5')

                ax.plot(coord_sdss_src.ra.value, coord_sdss_src.dec.value, transform=ax.get_transform('fk5'), marker="s",
                        markeredgecolor='orange', markerfacecolor='none', markersize=4, zorder=2)

        output = "{0:s}{1:s}_continuum_and_sdss.png".format(cfg_par['general'].get(
            'plotdir'), cfg_par['general']['workdir'].split('/')[-2])

        if cfg_par['general']['plot_format'] == "pdf":
            plt.savefig(output.replace(".png", ".pdf"),
                        overwrite=True, bbox_inches='tight')
        else:
            plt.savefig(output, overwrite=True, bbox_inches='tight', dpi=300)

        fits_hdulist.close()

    print("# Running SDSS source finder ... Done")


def match_sdss_to_radio(cfg_par):
    """This function matches the SDSS sources to the radio sources found in the APERTIF image

    """

    # Getting files and sources
    # +++++++++++++++++++++++++

    key = "sdss_match"

    workdir = cfg_par['general']['workdir']

    # sdss source file
    sdss_src_cat_file = "{0:s}abs/{1:s}_sdss_src.csv".format(
        workdir, workdir.split('/')[-2])

    sdss_src_cat = Table.read(sdss_src_cat_file, format='csv')

    # radio source cat
    radio_src_cat_file = "{0:s}abs/mir_src_sharpener.csv".format(
        workdir, workdir.split('/')[-2])

    radio_src_cat = Table.read(radio_src_cat_file, format='csv')

    # output file
    radio_sdss_src_cat_file = "{0:s}abs/{1:s}_radio_sdss_src.csv".format(
        workdir, workdir.split('/')[-2])

    # Match radio and SDSS by going through the radio sources
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # number of radio sources
    n_radio_src = np.size(radio_src_cat['ra'])

    # add columns
    radio_src_cat['sdss_id'] = Column(np.zeros(n_radio_src))
    radio_src_cat['sdss_ra'] = Column(np.zeros(n_radio_src))
    radio_src_cat['sdss_dec'] = Column(np.zeros(n_radio_src))
    radio_src_cat['sdss_radio_sep'] = Column(np.zeros(n_radio_src))
    radio_src_cat['sdss_redshift'] = Column(np.zeros(n_radio_src))

    # create coordinate objects
    coord_radio_src = SkyCoord(radio_src_cat['ra'], radio_src_cat['dec'], unit=(
        u.hourangle, u.deg), frame='fk5')
    coord_sdss_src = SkyCoord(
        sdss_src_cat['ra'], sdss_src_cat['dec'], unit=(u.deg, u.deg), frame='fk5')

    print("Matching radio sources")

    for k in range(n_radio_src):

        if radio_src_cat[k]['flux_int'] > cfg_par[key]['min_radio_flux']:

            idx, d2d, d3d = coord_radio_src[k].match_to_catalog_sky(
                coord_sdss_src)

            if d2d < cfg_par[key]['max_sep'] * u.arcsec:

                radio_src_cat[k]['sdss_id'] = sdss_src_cat[idx]['objid']
                radio_src_cat[k]['sdss_ra'] = coord_sdss_src[idx].ra.value
                radio_src_cat[k]['sdss_dec'] = coord_sdss_src[idx].dec.value
                radio_src_cat[k]['sdss_radio_sep'] = d2d[0].deg
                radio_src_cat[k]['sdss_redshift'] = sdss_src_cat[idx]['z']

    # save file
    radio_src_cat.write(radio_sdss_src_cat_file, format="csv", overwrite=True)
