__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys,string,os,math
import numpy as np
from astropy import units
from astropy.io import fits
from astropy import wcs
import logging


#define constants
RAD2DEG=180./math.pi
HI=1.42040575177e+09 #Hz
TSPIN=100            #K
MSUN=1.98855e33      #g
MHI=1.6749E-24       #g
CHI=2.36E5
PC=3.08567758E18    #cm
JANSKY=1e-23        #erg/scm2Hz
C=2.99792458E10     #cm/s
G=6.6742E-08        #cm3kg-1s-1      
MP=1.67492728E-24   #g
SIGMAT=6.66524E-25  #cm2

def ra2deg(ra_hms):
	'''
	Converts right ascension in hms coordinates to degrees and radians

	INPUT
	
	rahms: ra in HH:MM:SS format (str)
	
	OUTPUT
	
	conv_units.radeg: ra in degrees
	conv_units.rarad: ra in radians
	'''

	ra = string.split(ra_hms, ':')

	hh = float(ra[0])*15
	mm = (float(ra[1])/60)*15
	ss = (float(ra[2])/3600)*15

	return hh+mm+ss

def ra2hms( ra):

	RA , rs = '', ''
	if str(ra)[0] == '-':
		rs, ra = '-', abs(ra)
	raH = int(ra/15)
	raM = int(((ra/15)-raH)*60)
	raS = ((((ra/15)-raH)*60)-raM)*60
	RA = '{0}{1}:{2}:{3}'.format(rs, raH, raM, raS)		
	
	return ra

def dec2deg(dec_dms):

	'''
	Converts right ascension in hms coordinates to degrees and radians
	
	INPUT
	
	rahms: ra in HH:MM:SS format (str)
	
	OUTPUT
	
	conv_units.radeg: ra in degrees
	conv_units.rarad: ra in radians
	'''

	dec = string.split(dec_dms, ':')

	hh = abs(float(dec[0]))
	mm = float(dec[1])/60
	ss = float(dec[2])/3600

	if float(dec[0])>= 0:
		return hh+mm+ss
	else:
		return -(hh+mm+ss)

	return hh+mm+ss   

def dec2dms(ra='', dec='', round=False):
	
	DEC, ds = '', ''
  	if str(dec)[0] == '-':
	  ds, dec = '-', abs(dec)
	deg = int(dec)
	decM = abs(int((dec-deg)*60))
	decS = (abs((dec-deg)*60)-decM)*60
	DEC = '{0}{1}:{2}:{3}'.format(ds, deg, decM, decS)

	return DEC

def fov_of_cube(cubename):
	'''
	
	Module called by abs_ex
	Converts ra,dec of continuum sources
	into pixel coordinates of the datacube
	
	'''

	#I load the WCS coordinate system:
	#open file

	hdulist = fits.open(cubename)  # read input
	
	# read data and header
	#what follows works for wcs, but can be written better
	prihdr = hdulist[0].header
	if prihdr['NAXIS'] == 4:
		del prihdr['CTYPE4']
		del prihdr['CDELT4']    
		del prihdr['CRVAL4']
		del prihdr['CRPIX4']
	del prihdr['CTYPE3']
	del prihdr['CDELT3']
	del prihdr['CRVAL3']
	del prihdr['CRPIX3'] 
	del prihdr['NAXIS3']
	del prihdr['NAXIS']
	prihdr['NAXIS']=2
	

	w=wcs.WCS(prihdr) 

	foot = w.calc_footprint(prihdr)

	return foot

def coord_to_pix(imagename,ra,dec,verbose=False):
	'''
	
	Module called by abs_ex
	Converts ra,dec of continuum sources
	into pixel coordinates of the datacube
	
	'''

	#I load the WCS coordinate system:
	#open file

	hdulist = fits.open(imagename)  # read input

	# read data and header
	#what follows works for wcs, but can be written better
	# RS added some additional if clauses
	prihdr = hdulist[0].header
	if prihdr['NAXIS'] == 4:
		if 'CTYPE4' in prihdr:
			del prihdr['CTYPE4']
		if 'CDELT4' in prihdr:
			del prihdr['CDELT4']
		if 'CRVAL4' in prihdr:
			del prihdr['CRVAL4']
		if 'CRPIX4' in prihdr:
			del prihdr['CRPIX4']
		if 'CUNIT4' in prihdr:
			del prihdr['CUNIT4']
	
	if 'CTYPE3' in prihdr:
		del prihdr['CTYPE3']
	if 'CDELT3' in prihdr:
		del prihdr['CDELT3']
	if 'CRVAL3' in prihdr:
		del prihdr['CRVAL3']
	if 'CRPIX3' in prihdr:
		del prihdr['CRPIX3'] 
	if 'NAXIS3' in prihdr:
		del prihdr['NAXIS3']
	if 'CUNIT3' in prihdr:
		del prihdr['CUNIT3']

	del prihdr['NAXIS']
	prihdr['NAXIS']=2
	
	w=wcs.WCS(prihdr)    

	pixels=np.zeros([len(ra),2])
	count_out = 0
	count_flag = 0 
	for i in xrange(0,len(ra)):
		if ra[i] == 'nan':
			pixels[i, 0]= np.nan
			pixels[i, 1]= np.nan
			count_flag +=1
			if verbose == True:
				print('# Source # '+str([i])+ ' is flagged #')
		else:
			ra_deg = ra2deg(ra[i])
			dec_deg = dec2deg(dec[i])
			px,py=w.wcs_world2pix(ra_deg,dec_deg,0)
			if (0 < round(px,0) < prihdr['NAXIS1'] and
					0 < round(py,0) < prihdr['NAXIS2']): 
				pixels[i, 0]= round(px,0)
				pixels[i, 1]= round(py,0)
			else :
				pixels[i, 0]= np.nan
				pixels[i, 1]= np.nan
				count_out +=1
				if verbose == True:
					print '# Source # '+str([i])+ ' lies outside the fov of the data cube #'

	print '# Total number of sources: \t'+str(len(ra))
	print '# Sources below threshold: \t'+str(count_flag)
	print '# Sources outside f.o.v.:\t'+str(count_out)
	print '# Sources to analyze: \t\t'+str(len(ra)-count_flag-count_out)+'\n'

	return pixels

