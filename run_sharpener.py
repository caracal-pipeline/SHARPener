#!/usr/bin/python2

"""Script to run SHARPener

TODO:
1. Correct the name of the output plot
2. make script parallel

"""

import os
import string
import sys
import numpy as np
import pyaml
import json
import argparse
from astropy.io import ascii
import sharpener as sharpy
import tabulate
import cont_src as cont_src
import spec_ex as spec_ex
import time
import absorption_plot as abs_pl
import glob
from PyPDF2 import PdfFileMerger
import zipfile

# starting time of script
time_start = time.time()

# Create and parse argument list
# ++++++++++++++++++++++++++++++
parser = argparse.ArgumentParser(
    description='Run sharpener script')

# 1st argument: File name
parser.add_argument("yml_name", help='The setting file')

# 2nd-4th argument are optional: Enable or Disable certain steps
parser.add_argument("--do_cont", action="store_true", default=False,
                    help='Enable continuum extraction')

parser.add_argument("--do_spectra", action="store_true", default=False,
                    help='Enable spectra extraction')

parser.add_argument("--do_plots", action="store_true", default=False,
                    help='Enable plotting')

args = parser.parse_args()


# Set script parameters
# +++++++++++++++++++++
do_continuum_extraction = args.do_cont
do_spectra_extraction = args.do_spectra
do_plots = args.do_plots

# Load parameter file
# +++++++++++++++++++
time_start_load = time.time()

print("## Load parameter file")

spar = sharpy.sharpener(
    args.yml_name)

print("## Load parameter file ... Done ({0:.2f}s)".format(
    time.time() - time_start_load))

# print(spar.cfg_par)

# Find continuum sources
# ++++++++++++++++++++++
if do_continuum_extraction:
    time_start_find = time.time()

    print("## Find continuum sources")

    # get sources in continuum image
    sources = cont_src.find_src_imsad(spar.cfg_par)

    print("## Find continuum sources ... Done ({0:.2f}s)".format(
        time.time() - time_start_find))

    # print results
    src_list = ascii.read(spar.cfg_par['general']
                          ['absdir']+'mir_src_sharpener.csv')
    print(src_list)

# Extract spectra
# +++++++++++++++
if do_spectra_extraction:
    time_start_extract = time.time()

    print("## Extract HI spectra from cube")

    spectra = spec_ex.abs_ex(spar.cfg_par)

    print("## Extract HI spectra from cube ... Done ({0:.2f}s)".format(
        time.time() - time_start_extract))

# Plot spectra
# ++++++++++++
if do_plots:
    time_start_plot = time.time()

    print("## Plotting spectra")

    abs_pl.create_all_abs_plots(spar.cfg_par)

    print("## Plotting spectra ... Done ({0:.2f}s)".format(
        time.time() - time_start_plot))

# Merge plots:
# ++++++++++++
if spar.cfg_par['general']['merge_plots'] and spar.cfg_par['general']['plot_format'] == "pdf":

    time_start_merge = time.time()

    print("## Merging plots")

    # Merge the detailed plots
    # ++++++++++++++++++++++++
    plot_list = glob.glob(
        "{0:s}/*J*_detailed.pdf".format(spar.cfg_par['general']['plotdir']))

    if len(plot_list) != 0:

        plot_list.sort()

        plot_list.insert(0, "{0:s}{1:s}_continuum.pdf".format(
            spar.cfg_par['general']['plotdir'], spar.cfg_par['general']['workdir'].split("/")[-2]))

        pdf_merger = PdfFileMerger()

        for files in plot_list:
            pdf_merger.append(files)

        plot_name = "{0:s}{1:s}_all_plots_detailed.pdf".format(
            spar.cfg_par['general']['plotdir'], spar.cfg_par['general']['workdir'].split("/")[-2])

        pdf_merger.write(plot_name)
    else:
        print("No detailed plots found. Continue")

    # Merge compact plots
    # +++++++++++++++++++
    plot_list = glob.glob(
        "{0:s}/*J*_compact.pdf".format(spar.cfg_par['general']['plotdir']))

    if len(plot_list) != 0:

        plot_list.sort()

        plot_list.insert(0, "{0:s}{1:s}_continuum.pdf".format(
            spar.cfg_par['general']['plotdir'], spar.cfg_par['general']['workdir'].split("/")[-2]))

        pdf_merger = PdfFileMerger()

        for files in plot_list:
            pdf_merger.append(files)

        plot_name = "{0:s}{1:s}_all_plots_compact.pdf".format(
            spar.cfg_par['general']['plotdir'], spar.cfg_par['general']['workdir'].split("/")[-2])

        pdf_merger.write(plot_name)
    else:
        print("No plots found. Continue")

    print("## Merging plots ... Done ({0:.2f}s)".format(
        time.time() - time_start_merge))

    # Create a zip file with all the plots
    # ++++++++++++++++++++++++++++++++++++
    time_start_zip = time.time()
    print("## Create zip files")

    plot_list = glob.glob(
        "{0:s}/*all_plots*.pdf".format(spar.cfg_par['general']['plotdir']))

    if len(plot_list) != 0:

        with zipfile.ZipFile('NGC315_plots.zip', 'w') as myzip:

            for plot in plot_list:
                myzip.write(plot, os.path.basename(plot))
    else:
        print("No files to zip.")

    print("## Create zip files ... Done ({0:.2f}s)".format(
        time.time() - time_start_zip))
# Script finished
print("Script finished ({0:.2f}s)".format(
    time.time() - time_start))