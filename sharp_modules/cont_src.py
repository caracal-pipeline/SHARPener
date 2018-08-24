__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, string, os
import numpy as np
import yaml, glob
import convert_units as conv_units
from astropy import wcs
from astropy.io import fits as pyfits
from astropy.io import ascii
from astropy import units as u
from astropy.table import Table, Column, MaskedColumn
from astroquery.vizier import Vizier
import astropy.coordinates as coord



import mirlib as lib
import convert_units as conv_units


####################################################################################################



def source_catalog(cfg_par,tablename):
	'''
	Takes a dictionary with the following parameters:


	'''
	key = 'source_catalog'
	width = cfg_par[key].get('width', '2')
	catalog = cfg_par[key].get('catalog', 'NVSS')
	thresh = cfg_par[key].get('thresh', 10e-3)
	centre = cfg_par[key].get('centre_coord',None)
	catalog_table = tablename 
	Vizier.ROW_LIMIT = -1

	p = Vizier.query_region(coord.SkyCoord(centre[0],centre[1], unit=(u.hourangle, u.deg), 
		frame = 'icrs'), width = width, catalog = catalog)
	tab = p[0]
	ra_deg = []
	dec_deg = []
	if catalog == 'NVSS':
		for i in xrange (0, len(tab['RAJ2000'])):
		   tab['RAJ2000'][i] = string.join(string.split(tab['RAJ2000'][i],' '),':')
		   ra_deg.append(conv_units.ra2deg(tab['RAJ2000'][i]))
		   tab['DEJ2000'][i] = string.join(string.split(tab['DEJ2000'][i],' '),':')
		   dec_deg.append(conv_units.dec2deg(tab['DEJ2000'][i]))

		above_thresh = tab['S1.4']<thresh
	for i in xrange(1,len(tab.colnames)):
		tab[tab.colnames[i]][above_thresh] = np.nan

	tab =  Table(tab, masked=True)
	ascii.write(tab, catalog_table, overwrite=True)

	return tab

def sim_cont_from_cube(tablename,catalog,infile,outfile):

	tab = ascii.read(table)

	if catalog == 'NVSS':

			ra = tab['RAJ2000']
			dec = tab['DEJ2000']
	
			major = tab['MajAxis']
			minor = tab['MinAxis']		

			angle1=np.radians(0.0)
			cosangle1=np.cos(angle1)
			sinangle1=np.sin(angle1)

	pixels = conv_units.coord_to_pix(infile,ra,dec,verbose=False)


	cubefile = pyfits.open(infile)
	contdata = cubefile[0].data
	cubehead = cubefile[0].header

	if cubehead['NAXIS'] > 3:

		contdata = np.zeros([contdata.shape[2],contdata.shape[3]])
		if 'CTYPE4' in cubehead:
			del cubehead['CTYPE4']
		if 'CDELT4' in cubehead:
			del cubehead['CDELT4']    
		if 'CRVAL4' in cubehead:
			del cubehead['CRVAL4']
		if 'CRPIX4' in cubehead:
			del cubehead['CRPIX4']
		if 'NAXIS4' in cubehead:
			del cubehead['NAXIS4']

	elif cubehead['NAXIS'] == 3:

		contdata = np.zeros([contdata.shape[1],contdata.shape[2]])
		if 'CRPIX3' in cubehead:
			del cubehead['CRPIX3'] 
		if 'CRVAL3' in cubehead:
			del cubehead['CRVAL3']
		if 'CDELT3' in cubehead:
			del cubehead['CDELT3']
		if 'CTYPE3' in cubehead:
			del cubehead['CTYPE3']
		if 'NAXIS3' in cubehead:
			del cubehead['NAXIS3']

	else: 
		contdata = np.zeros([contdata.shape[0],contdata.shape[1]])



	del cubehead['NAXIS']

	w=wcs.WCS(cubehead)    

	xnum = np.linspace(0,cubehead['NAXIS1'],cubehead['NAXIS1'])
	ynum = np.linspace(0,cubehead['NAXIS2'],cubehead['NAXIS2'])	
	x, y = 	np.meshgrid(xnum, ynum)

	for i in xrange(0,pixels.shape[0]):

		xc=float(pixels[i][0])
		yc=float(pixels[i][1])

		if str(xc) == 'nan' or str(yc) == 'nan':
			pass
		else:
			if minor[i]/3600. >= float(cubehead['CDELT2']) and major[i]/3600. >= float(cubehead['CDELT2']):
				a = major[i]/3600./float(cubehead['CDELT2'])/2.
				b = minor[i]/3600./float(cubehead['CDELT2'])/2.

				ell = np.power(x-xc, 2)/np.power(a,2) + np.power(y-yc, 2)/np.power(b,2)
				index_ell = np.where(np.less_equal(ell,1))
				contdata[index_ell] = 1
			else:
				contdata[int(yc),int(xc)] = 1

	pyfits.writeto(outfile,contdata,cubehead,overwrite=True)

def write_src_csv(tot,cfg_par):
	'''
		- This module saves the output of MIRIAD in a csv table
	INPUT:
		tot: array with output of imsad (final product of find_src_imsad)
	OUTPUT
		mir_src_sharpener.csv : default name, table saved in absdir (see parameter file)
	'''

	# write the spectrum on file
	out_file = str(cfg_par['general'].get('absdir')) + 'mir_src_sharpener.csv'

	f = open(out_file, 'w')
	f.write('ID,J2000,ra,dec,Pix_x,Pix_y,peak\n')
	np.savetxt(f, tot, delimiter=",", fmt="%s")
	f.close()
	
	print '# List of continuum sources saved on file. #' 

	return 0

def mosaic_continuum(cfg_par):
		
		os.chdir(cfg_par['general']['workdir'])

		print '# Convert continuum images to miriad format'   
	
		key = 'mosaic_continuum'
		mosaicdir = cfg_par[key]['mosaic_directory']


		#baseimage = cfg_par[key]['base_image']
		#baseimage_mir = string.split(baseimage,'.')[0]
		#baseimage_mir = baseimage_mir+'.mir'
		#mirfits = lib.miriad('fits')
		#mirfits.op = 'xyin'
		#mirfits.in_ = baseimage
		#mirfits.out = baseimage_mir
		#mirfits.go(rmfiles=True)
		
		os.chdir(mosaicdir)
		filelist = [f for f in glob.glob('*.fits')]
		mirlist=[]
		for i in xrange(0,len(filelist)):		
			ff = pyfits.open(filelist[i])
			heads = ff[0].header
			ra_deg = float(heads['CRVAL1'])
			dec_deg = float(heads['CRVAL2'])
			rad = np.pi/180.*ra_deg
			decd = np.pi/180.*dec_deg

			mircont = os.path.basename(filelist[i])
			mircont = string.split(mircont,'.')[0]
			mircont = str(i)+'.mir'
			
			mirfits = lib.miriad('fits')
			mirfits.op = 'xyin'
			mirfits.in_ = filelist[i]
			mirfits.out = mircont
			mirfits.go(rmfiles=True)
			mirlist.append(mircont)

			puthd = lib.miriad('puthd')
			puthd.in_=mircont+'/pbfwhm'
			puthd.value=3409.5199738
			puthd.go()
			puthd = lib.miriad('puthd')

			puthd.type = 'double'
			puthd.in_=mircont+'/pntra'
			puthd.value=ra_deg
			puthd.type = 'r'
			puthd.go()
			puthd = lib.miriad('puthd')

			puthd.in_=mircont+'/pntdec'
			puthd.value=dec_deg
			puthd.type = 'r'
			puthd.go()
			puthd = lib.miriad('puthd')

			puthd.in_=mircont+'/pbinit'
			puthd.value=1.370031054688e9
			puthd.type = 'double'
			puthd.go()

		# print '# Regrid continuum images to base image'   
		# mireg_list=[]
		# for i in xrange(0,len(mirlist)):

		# 	mireg = os.path.basename(mirlist[i])
		# 	mireg = string.split(mireg,'.')[0]
		# 	mireg = str(i)+'_reg.mir'			

		# 	miregrid = lib.miriad('regrid')
		# 	miregrid.tin = '../'+baseimage_mir
		# 	miregrid.in_ = mirlist[i]
		# 	miregrid.bw = cfg_par[key].get('mosaic_bw',0)
		# 	miregrid.out = mireg
		# 	miregrid.go(rmfiles=True)
		# 	mireg_list.append(mireg)
		# 	puthd = lib.miriad('puthd')
		# 	puthd.in_=mireg+'/pbfwhm'
		# 	puthd.value=3409.5199738
		# 	puthd.go()
		# 	puthd = lib.miriad('puthd')

		# 	puthd.type = 'double'
		# 	puthd.in_=mireg+'/pntra'
		# 	puthd.value=ra_deg
		# 	puthd.type = 'r'
		# 	puthd.go()
		# 	puthd = lib.miriad('puthd')

		# 	puthd.in_=mireg+'/pntdec'
		# 	puthd.value=dec_deg
		# 	puthd.type = 'r'
		# 	puthd.go()
		# 	puthd = lib.miriad('puthd')

		# 	puthd.in_=mireg+'/pbinit'
		# 	puthd.value=1.370031054688e9
		# 	puthd.type = 'double'
		# 	puthd.go()
			
		

		# print '# Mosaic continuum images'   

		# mirmos = lib.miriad('linmos')
		# mirmos.tin = '../'+baseimage_mir
		# print mirmos.tin
		# mirmos.in_ = mireg_list
		# print mirmos.in_
		# mirmos.out = 'mosaic.mir'
		# mirmos.rms = cfg_par[key]['mosaic_rms']
		# print mirmos.rms
		# mirmos.cutoff = cfg_par[key]['mosaic_cutoff']
		# mirmos.go(rmfiles=True)	

		# fits = lib.miriad('fits')
		# fits.op = 'xyout'
		# fits.in_ = 'mosaic.mir'
		# fits.out = 'mosaic.fits'
		# fits.go(rmfiles=True)


def find_src_imsad(cfg_par):
	'''

	Finds the continuum sources according to the options set in the source_finder sub-keys
	INPUT:
		cfg_par : parameter file loaded with sharpener
	OUTPUT:
		csv table of continuum sources (output of MIRIAD imsad)

	'''		
	print '# Find continuum sources '   

	os.chdir(cfg_par['general']['workdir'])
	cont_im_mir = cfg_par['general']['mircontname']
	cont_im = os.path.basename(cfg_par['general']['contname'])

	src_imsad_out = 'abs/mir_src_sharpener.txt'
	key = 'source_finder'

	if os.path.exists(cont_im_mir) == False and os.path.exists(cont_im) == True: 
	
		fits = lib.miriad('fits')
		fits.op = 'xyin'
		fits.in_ = cont_im
		fits.out = cont_im_mir
		fits.go(rmfiles=True)

	elif os.path.exists(cont_im) == False and os.path.exists(cont_im_mir) == True: 
	
		fits = lib.miriad('fits')
		fits.op = 'xyout'
		fits.in_ = cont_im_mir
		fits.out = cont_im
		fits.go(rmfiles=True)

	
	imsad = lib.miriad('imsad')
	imsad.in_ =cont_im_mir
	imsad.out = src_imsad_out
	imsad.clip = cfg_par[key]['clip']
	#imsad.region = 'boxes\('+self.abs_ex_imsad_region+'\)'
	imsad.options = cfg_par[key]['options']

	imsad.go(rmfiles=True)
	
	
	#modify output of imsad for module load_src_csv
	src_list = open(src_imsad_out,'r')
	lines = src_list.readlines()
	len_lines = len(lines)
	
	ra_tmp = []
	dec_tmp = []
	peak_tmp = []        
	
	for i in xrange (0,len_lines):
		lines[i] = lines[i].strip()
		tmp = lines[i].split(' ')
		ra_tmp.append(str(tmp[3]))
		dec_tmp.append(str(tmp[4]))
		peak_tmp.append(str(tmp[6]))
	
	ra_tmp = np.array(ra_tmp)
	dec_tmp = np.array(dec_tmp)
	peak_tmp = np.array(peak_tmp)
   
	#ID
	ids = np.array(np.arange(1,len_lines+1,1),dtype=str)
	
	#J2000
	#convert ra
	ra_vec = []
	for i in xrange (0, len_lines):
		line = ra_tmp[i].split(':')
		last_dig = int(round(float(line[2]),0))
		if last_dig < 10:
			last_dig = '0'+str(last_dig)
		ra_vec.append(line[0]+line[1]+str(last_dig))
	
	#convert dec 
	dec_vec = []
	dec_coord = []
	for i in xrange (0, len_lines):
		line = dec_tmp[i].split(':')
		last_dig = int(round(float(line[2]),0))
		first_dig = line[0].split('-')
		dec_vec.append(first_dig[1]+line[1]+str(last_dig))
		dec_coord.append('-'+first_dig[1]+':'+line[1]+':'+str(last_dig))        

	J2000_tmp = np.array([ a+'+'+b for a,b in zip(ra_vec,dec_vec)])
	dec_coord = np.array(dec_coord)
   
	### Find pixels of sources in continuum image
	#open continuum image
	hdulist = pyfits.open(cont_im)  # read input
	# read data and header
	#what follows works for wcs, but can be written better
	prihdr = hdulist[0].header   
	if 'CTYPE4' in prihdr:
		del prihdr['CTYPE4']
	if 'CDELT4' in prihdr:
		del prihdr['CDELT4']    
	if 'CRVAL4' in prihdr:
		del prihdr['CRVAL4']
	if 'CUNIT4' in prihdr:
		del prihdr['CUNIT4']
	if 'CRPIX4' in prihdr:			
		del prihdr['CRPIX4']
	if 'CTYPE3' in prihdr:
		del prihdr['CTYPE3']
	if 'CDELT3' in prihdr:
		del prihdr['CDELT3']
	if 'CRVAL3' in prihdr:
		del prihdr['CRVAL3']
	if 'CRPIX3' in prihdr:
		del prihdr['CRPIX3'] 
	if 'CUNIT3' in prihdr:
		del prihdr['CUNIT3']
	if 'NAXIS3' in prihdr:
		del prihdr['NAXIS3']
	if 'NAXIS' in prihdr:
		del prihdr['NAXIS']
	prihdr['NAXIS']=2
	#load wcs system
	w=wcs.WCS(prihdr)    

	pixels_cont=np.zeros([ra_tmp.shape[0],2])
	ra_list_tmp = np.zeros([ra_tmp.shape[0]])
	for i in xrange (0,ra_tmp.shape[0]):
		
		ra_deg_tmp = conv_units.ra2deg(ra_tmp[i])
		dec_deg_tmp = conv_units.dec2deg(dec_coord[i])
		
		px,py=w.wcs_world2pix(ra_deg_tmp,dec_deg_tmp,0)
		pixels_cont[i,0]= str(round(px,0))
		pixels_cont[i,1]= str(round(py,0))        
	
	#make csv file with ID, J2000, Ra, Dec, Pix_y, Pix_y, Peak[Jy] 
	tot = np.column_stack((ids,J2000_tmp,ra_tmp,dec_coord,
						   pixels_cont[:,0],pixels_cont[:,1],peak_tmp))

	write_src_csv(tot,cfg_par)  
	
	print '# Continuum sources found #'    

	
	return tot  






