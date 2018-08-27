__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, string, os
import numpy as np
import yaml
#from continuum import *
import json
import convert_units as conv_units
import cont_src as cont_src
import cubez as cubef
from kk import *
import hi 

#import radiobs
#from radiobs import conv_units, cubeful, hi
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

hi = hi.hi()

C=2.99792458e5 #km/s
HI=1.420405751e9 #Hz

####################################################################################################

def abs_ex(cfg_par):

		verb = cfg_par['general']['verbose']

		cubename = cfg_par['general'].get('cubename',None)
		cubefile = fits.open(cubename)  # read input
		hdr = cubefile[0].header
		sci = cubefile[0].data 
		sci = sci.squeeze()
		x = hdr['NAXIS1']
		y = hdr['NAXIS2']
		z = hdr['NAXIS3']
		cen_imx = hdr['CRPIX1']
		cen_imy = hdr['CRPIX2']
		freq0 = hdr['CRVAL3']
		freq_del = hdr['CDELT3']

		key = 'source_catalog'
		if cfg_par['source_catalog'].get('enable',False) == True:

			catalog_table = str(cfg_par['general'].get('absdir')) + 'cat_src_sharpener.txt'
		   
			tab = ascii.read(catalog_table)
		
			if cfg_par[key].get('catalog', 'NVSS') == 'NVSS':

				J2000_name = tab['NVSS']
				ra = tab['RAJ2000']
				dec = tab['DEJ2000']
				flux_cont = tab['S1.4']*1e-3

		elif cfg_par['source_finder'].get('enable',False) == True:

			src_list_csv = cfg_par['general']['absdir']+'mir_src_sharpener.csv'
			# open file
			src_list_vec = ascii.read(src_list_csv)
			J2000_name = np.array(src_list_vec['J2000'],dtype=str)
			ra = np.array(src_list_vec['ra'],dtype=str)
			dec = np.array(src_list_vec['dec'],dtype=str)
			flux_cont = np.array(src_list_vec['peak'],dtype=float)
		


		pixels = conv_units.coord_to_pix(cubename,ra,dec, verbose=False)

		key = 'spec_ex'

		src_id = np.arange(0,ra.size+1,1)
		freq = cubef.zaxis(cubename)

		abs_mean_rms = np.zeros(pixels.shape[0])
		abs_los_rms = np.zeros(pixels.shape[0])
		tau_los_rms = np.zeros(pixels.shape[0])
		outnames = []
		count_thresh =0
		count_fov = 0
		count_blanks = 0
		average_noise = []
		for i in xrange(0,pixels.shape[0]):

			# extract spectrum from each line of sight
			flux = np.zeros(freq.shape[0])
			madfm = np.zeros(freq.shape[0])

			if str(pixels[i,0]) == 'nan' or str(pixels[i,1]) == 'nan':
				count_thresh +=1
				pass

			elif (0 < int(pixels[i,0]) < x and
					0 < int(pixels[i,1]) < y): 
					
				pix_x_or = int(pixels[i,0])
				pix_y_or = int(pixels[i,1])
				for j in xrange(0, z):
					chrom_aber = cfg_par[key].get('chrom_aberration', False)
					#correct for chromatic aberration
					if chrom_aber == True:

						if (cfg_par[key].get('zunit','Hz') == 'm/s'):
							freq_real= freq* 1e2
							freq_real = (kk.C*kk.HI) /  (freq + kk.C)
							freq_real0 = (kk.C*kk.HI) /  (hdr['CRVAL3']*1e2 + kk.C)
							freq_del = (freq_real0 - freq_real[-1] )/ len(freq_real)
						#depending if the cube is in velocity or frequency ?
							scale = (freq_real0 - j*freq_del) / freq_real0

						pix_x = (pix_x_or - hdr['CRPIX1']) * scale + hdr['CRPIX1']
						pix_y = (pix_y_or - hdr['CRPIX1']) * scale + hdr['CRPIX1']
						pix_x = int(round(pix_x,0))
						pix_y = int(round(pix_y,0))
					else:
						pix_x = pix_x_or
						pix_y = pix_y_or
					
					if  (0 < pix_x < x and
						 0 < pix_y < y): 
						flux[j] = sci[j, pix_y, pix_x]
					else:
						flux[j] = 0.0


					# determine the noise of the spectrum [Whiting 2012 et al.] in each channel
					# MADMF: median absolute deviation from the median
					# extract a region were to determine the noise: A BOX around the l.o.s.
					if (pix_x+10 < hdr['NAXIS1'] and  pix_x-10 > 0 and
					   pix_y+10 < hdr['NAXIS2'] and pix_y - 10 > 0):
							rms = np.nanmedian(sci[j, pix_x +10:pix_x + 10, pix_y - 5:pix_y + 5])
							if rms != 0.0:
								med2 = np.abs(sci[j, pix_y, pix_x] - rms)
								madfm[j] = np.nanmedian(med2) / 0.6744888
							else:
								madfm[j] = 0.0
					else:
						madfm[j] = 0.0


					abs_mean_rms[i] = np.nanmean(madfm) 

				if np.nansum(flux) == 0.:
					count_blanks +=1
					if verbose == True:
						print '# Blank spectrum:\t'+str(src_id[i])+' '+J2000_name[i]+' #'
					continue

				# measure noise in the spectrum outside of the line
				end_spec = float(sci.shape[0])
				end_spec_th = int(end_spec/3.)
				end_spec = int(end_spec)
				mean_rms = (np.std(flux[0:end_spec_th]) +
									np.std(flux[end_spec-end_spec_th:end_spec])) / 2.
				mean_rms_arr = np.zeros(sci.shape[0])+mean_rms
				
				average_noise.append(mean_rms)

				tau = hi.optical_depth(flux, flux_cont[i])
				if np.nansum(madfm)!= 0.0:
					tau_noise = hi.optical_depth(madfm, flux_cont[i])
				else:
					tau_noise = np.zeros(sci.shape[0])

				#write spectrum
				out_spec = str(cfg_par['general']['specdir']+str(src_id[i])+'_J'+J2000_name[i])+'.txt'
				outnames.append(out_spec)

				flag_chans = cfg_par[key].get('flag_chans', None)
				if flag_chans != None:
					index_flags_l = (np.abs(freq - flag_chans[0])).argmin()
					for k in xrange(1,len(flag_chans)):
						index_flags = (np.abs(freq - flag_chans[k])).argmin()
						flux[index_flags_l:index_flags] = np.nan
						index_flags_l = index_flags

				if cfg_par[key].get('zunit','Hz') == 'm/s':
					xcol = 'Velocity [m/s]'
				elif cfg_par[key].get('zunit','Hz') == 'km/s':
					xcol = 'Velocity [km/s]'
				elif cfg_par[key].get('zunit','Hz') == 'MHz':
					xcol = 'Frequency [MHz]'
				else:
					xcol = 'Frequency [Hz]'

				t = Table([freq, flux, madfm, tau, tau_noise, mean_rms_arr], 
					names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
					meta={'name': 'Spectrum'})
				ascii.write(t,out_spec,overwrite=True)
				if verb==True:
					print '# Extracted spectrum: \t' +str(src_id[i])+' '+J2000_name[i]+' #'

				polysub = cfg_par['polynomial_subtraction'].get('enable', False) 
				if polysub == True:

					deg = cfg_par['polynomial_subtraction'].get('degree_pol',3)
					sub_flux = poly_sub(cfg_par,freq,flux,deg)
					sub_madfm = madfm.copy()
					sub_od = tau.copy()
					sub_noise_od = tau_noise.copy()

					out_spec_polysub = string.split(out_spec,'.')[0]
					out_spec_polysub= out_spec_polysub+'_psub.txt'

					t = Table([freq, sub_flux, sub_madfm, sub_od, sub_noise_od, mean_rms_arr], 
						names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
						meta={'name': 'Spectrum'})
					ascii.write(t,out_spec_polysub,overwrite=True)

				dohan = cfg_par['hanning'].get('enable', False)
				if dohan == True and polysub == False:
					
					window = cfg_par['hanning'].get('window', 1)
					han_flux  = hanning_spec(flux)
					han_madfm = hanning_spec(madfm)
					han_od = hanning_spec(tau)
					han_noise_od = hanning_spec(tau_noise)
					
					out_spec_han = string.split(out_spec,'.')[0]
					out_spec_han= out_spec_han+'_han.txt'

					t = Table([freq, han_flux, han_madfm, han_od, han_noise_od, mean_rms_arr], 
						names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
						meta={'name': 'Spectrum'})
					ascii.write(t,out_spec_han,overwrite=True)

				elif polysub == True and dohan == True:

					han_sub_flux  = hanning_spec(sub_flux)
					han_sub_madfm = hanning_spec(sub_madfm)
					han_sub_od = hanning_spec(sub_od)
					han_sub_noise_od = hanning_spec(sub_noise_od)
					
					out_spec_han = string.split(out_spec,'.')[0]
					out_spec_han= out_spec_han+'psub_han.txt'


					t = Table([freq, han_sub_flux, han_sub_madfm, han_sub_od, han_sub_noise_od, mean_rms_arr], 
						names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
						meta={'name': 'Spectrum'})
					ascii.write(t,out_spec_han,overwrite=True)




		print '\n# Total number of sources: \t'+str(pixels.shape[0])
		print '# Sources flagged: \t\t'+str(count_thresh)
		print '# Blank spectra:\t\t'+str(count_blanks)
		print '# Total number of spectra: \t'+str(pixels.shape[0]-count_thresh-count_fov-count_blanks)
		print '# Average noise in spectra: \t'+str(round(np.nanmean(average_noise)*1e3,1))+' mJy/beam'

		return 0 

def hanning_spec(flux):

		new_flux = flux.copy()
		new_flux[0] = (flux[0]+flux[1])/2.
		for i in xrange(1,len(flux)-1):

			new_flux[i] = (flux[i-1]+2.*flux[i]+flux[i+1])/4.

		new_flux[-1] = (flux[-2]+flux[-1])/2.

		return new_flux

def poly_sub(cfg_par,x, y,deg):

		if cfg_par['abs_ex'].get('zunit') == 'm/s':
			unit_z = u.m / u.s
		else:
			unit_z = u.Hz

		step = (x[-1]-x[0])/len(x)
		wave_for_spec = WaveCoord(cdelt=step, crval=x[0], cunit= unit_z)
		spe = Spectrum(wave=wave_for_spec, data=y)
		cont = spe.poly_spec(deg)

		cont_sub = y - cont.data



		return cont_sub



