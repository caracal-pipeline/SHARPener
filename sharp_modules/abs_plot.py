__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, string, os
import numpy as np
import yaml
import json
from astropy import wcs
from astropy.io import fits, ascii
from astropy import units as u
from astropy.table import Table, Column, MaskedColumn
from astroquery.vizier import Vizier
import astropy.coordinates as coord
from mpdaf.obj import Spectrum, WaveCoord
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import rc


C=2.99792458e5 #km/s
HI=1.420405751e9 #Hz

####################################################################################################


def abs_plot(spec_name,cfg_par):
		'''		
		Plots spectra of all radio sources found by find_src_imsad 
		saved in basedir/beam/abs/spec.
		Plots are stored in basedir/beam/abs/plot
		
		IN
			Spectra extracted by spec_ex
		
		IN cfga
			abs_ex_plot_xaxis= ' '      #: X-axis units ['velocity','frequency'] 
			abs_ex_plot_yaxis= ' '      #: Y axis units ['flux','optical depth']
			abs_ex_plot_redsrc= True    #: plots line at redshift of source in spectrum redshift must be stored in table of load_src_csv
			abs_ex_plot_title= True     #: plot title: J2000 name of radio source
			abs_ex_plot_format= ' '     #: format of plot ['.pdf','.jpeg','.png']
		
		OUT
			For each source outputs have the following name:
			J2000_xaxis-unit_yaxis-unit.plot_format = J220919.87+180920.17_vel_flux.pdf

		'''

		key = 'abs_plot'

		os.chdir(cfg_par['general']['specdir'])
		
		params = {
				  'text.usetex': True,
				  'text.latex.unicode': True
				   }
		rc('font', **{'family': 'serif', 'serif': ['serif']})        
		plt.rcParams.update(params)
		
		#for i in xrange(0,len(np.atleast_1d(spec_src_name))):
			
			#load data and labels 
		#	spec_name = spec_src_name[i]
		#	print spec_name
		if os.path.isfile(spec_name) == True:
		
			# Set plot specs
			font_size = 16            
			plt.ioff()
			fig = plt.figure(figsize =(9,6))
			fig.subplots_adjust(hspace=0.0)
			gs = gridspec.GridSpec(1, 1)
			plt.rc('xtick', labelsize=font_size-2)
			plt.rc('ytick', labelsize=font_size-2) 

			# Initialize subplots
			ax1 = fig.add_subplot(gs[0])
			ax1.set_xlabel('')
			ax1.set_ylabel('') 
			
		
			spec_vec = ascii.read(spec_name)


			x_data = np.array(spec_vec[spec_vec.colnames[0]],dtype=float)

			flag_chans = cfg_par['abs_ex'].get('flag_chans', None)
			flag_chans1 = cfg_par['abs_ex'].get('flag_chans', None)

			if cfg_par['abs_ex'].get('zunit') == 'm/s':
				x_data /= 1e3
				ax1.set_xlabel(r'$cz\,(\mathrm{km}\,\mathrm{s}^{-1})$', fontsize=font_size)


			y_data = np.array(spec_vec[spec_vec.colnames[1]],dtype=float)*1e3
			y_sigma = np.array(spec_vec[spec_vec.colnames[2]])

			if cfg_par['abs_ex'].get('zunit') == 'MHz':
				x_data /= 1e6
				ax1.set_xlabel(r'Frequency [MHz]', fontsize=font_size)
			

			ylabh = ax1.set_ylabel(r'S\,$[\mathrm{mJy}\,\mathrm{beam}^{-1}]$', fontsize=font_size)
			ylabh.set_verticalalignment('center')


			# Plot spectra 
#                if self.abs_ex_plot_linestyle == 'step':
				#ax1.plot(x_data, y_data, color='black', linestyle='-')

			if flag_chans != None:
				flag_chans = np.array(flag_chans) 
				if 	cfg_par['abs_ex'].get('zunit') == 'm/s':
					flag_chans = np.divide(flag_chans,1e3)
				index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
				for k in xrange(1,len(flag_chans)):
					index_flags = (np.abs(x_data - flag_chans[k])).argmin()
					#y_data[index_flags] = 0.0
				y_data[index_flags_l:index_flags] = 0.0
			ax1.step(x_data, y_data, where='mid', color='black', linestyle='-')


			# Calculate axis limits and aspect ratio
			x_min = np.min(x_data)
			x_max = np.max(x_data)
			y1_array = y_data[np.where((x_data>x_min) & (x_data<x_max))]
			y1_min = np.nanmin(y_data)*1.1
			y1_max = np.nanmax(y_data)*1.1

			# Set axis limits
			ax1.set_xlim(x_min, x_max)
			ax1.set_ylim(y1_min, y1_max)
			ax1.xaxis.labelpad = 6
			ax1.yaxis.labelpad = 10 


			if flag_chans1 != None:
				flag_chans = np.array(flag_chans1) 
				if 	cfg_par['abs_ex'].get('zunit') == 'm/s':
					flag_chans = np.divide(flag_chans,1e3)
				index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
				for k in xrange(1,len(flag_chans)):
					index_flags = (np.abs(x_data - flag_chans[k])).argmin()
					ax1.fill_between([x_data[index_flags_l],x_data[index_flags]],y1_min,y1_max, facecolor='grey', alpha=0.3)


			# Plot noise
			ax1.fill_between(x_data, -y_sigma, y_sigma, facecolor='grey', alpha=0.5)




			# Plot stuff
			ax1.axhline(color='k', linestyle=':', zorder=0)
			
			redshifts = cfg_par[key].get('redshift_sources',None)
			if len(redshifts) == 2:
					ax1.fill_between([redshifts[0],redshifts[1]],y1_min,y1_max, facecolor='red', alpha=0.1)

			# Add title        
			#if self.abs_ex_plot_title == True:
			#	ax1.set_title('%s' % (self.J2000_name[i]), fontsize=font_size+2) 
			#ax1.axes.titlepad = 8

			# Add minor tick marks
			ax1.minorticks_on()

			# Save figure to file
			outplot = os.path.basename(spec_name)
			outplot = string.split(outplot,'.')[0]
			outplot = cfg_par['general']['plotdir']+outplot+'.png'
			plt.savefig(outplot,
						overwrite = True)
			print '# Plotted spectrum of source ' + os.path.basename(spec_name)+'. #'
		else:
			print '# Missing spectrum of source ' + os.path.basename(spec_name)+'. #'