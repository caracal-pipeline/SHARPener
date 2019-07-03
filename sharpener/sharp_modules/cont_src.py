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
#from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, MaskedColumn, hstack
from astroquery.vizier import Vizier
import astropy.coordinates as coord

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mc

import mirlib as lib
import convert_units as conv_units
import absorption_plot as abs_plot

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
    catalog_table = cfg_par['general']['workdir']+tablename
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

    tab = ascii.read(tablename)

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
    x, y =  np.meshgrid(xnum, ynum)

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

    return 0

def write_src_csv(tot,cfg_par):
    """This funcions saves the output of MIRIAD in a csv table

    Parameter
    ---------
    tot: array with output of imsad (final product of find_src_imsad)
    
    Return
    ------
    mir_src_sharpener.csv : default name, table saved in absdir (see parameter file)
    """

    # write the spectrum on file
    out_file = str(cfg_par['general'].get('absdir')) + 'mir_src_sharpener.csv'

    f = open(out_file, 'w')
    f.write('ID,J2000,ra,dec,Pix_x,Pix_y,peak,flux_int,beam_major_decon,beam_minor_decon,beam_pang_decon,FLAG,DFLAG,FFLAG\n')
    np.savetxt(f, tot, delimiter=",", fmt="%s")
    f.close()
    
    print('# List of continuum sources saved on file. #')

    return 0

def check_miriad_output(imsad_file):
    """Function to check the imsad output file

    The function corrects the lengths of the imsad
    output file from miriad assuming that it is basically
    the last flag and not more. If imsad does not produce
    an FFLAG entry for some sources (or all), then it is left
    without value. This creates different column numbers for
    different rows.

    Parameter:
        imsad_file : str
            Name of imsad output file

    Return:
        No return value
    
    Output:
        Overwrites existing imsad file with corrected version
    """

    print("Checking imsad output file from Miriad")

    with open(imsad_file, 'r') as f:
        src_list = f.readlines()

    # get the length of the lines
    line_length_list = np.array([len(line) for line in src_list])
    min_line_length = np.min(line_length_list)
    max_line_length = np.max(line_length_list)

    if min_line_length != max_line_length:
        print("Found lines with different lengths. Attempting to correct.")

        if max_line_length == 98 and min_line_length == 96:
            short_lines_index_list = np.where(line_length_list == 96)[0]
            for line_index in short_lines_index_list:
                src_list[line_index] = src_list[line_index].replace("\n", " -\n")
        else:
            print("ERROR: Unknown difference in columns of imsad output")
            sys.exit(1)

        # save the file temporarily
        curPath = os.getcwd()
        tmp_file = curPath +"/sharpOut/abs/tmp.txt"
        with open(tmp_file, "w") as f:
            f.writelines(src_list)

        # rename temporary file by overwriting original file
        os.rename(tmp_file, imsad_file)

        print("... Done. Continue")
    # in case there are only 96 lines
    elif min_line_length == max_line_length and max_line_length == 96:
        print("Found lines with different lengths. Attempting to correct.")

        short_lines_index_list = np.where(line_length_list == 96)[0]
        for line_index in short_lines_index_list:
            src_list[line_index] = src_list[line_index].replace("\n", " -\n")

        # save the file temporarily
        curPath = os.getcwd()
        tmp_file = curPath +"/sharpOut/abs/tmp.txt"
        with open(tmp_file, "w") as f:
            f.writelines(src_list)

        # rename temporary file by overwriting original file
        os.rename(tmp_file, imsad_file)

        print("... Done. Continue")
    else:
        print("... File looks good. Continue.")

def create_karma_annotation_file(coord_list,cfg_par):
    """Function to create carma annotation file of the extracted source

    The function uses the astropy table objects to create the karma
    annotation file
    """

    print("Creating karma annotation file")

    with open(cfg_par['general']['absdir']+'karma_src_sharpener.ann', "w") as f:
        f.write("COORD W\n")
        f.write("PA STANDARD\n")
        f.write("COLOR GREEN\n")
        for coords in coord_list:
            f.write("CIRCLE {0:.5f} {1:.5f} {2:.5f}\n".format(
                coords.ra.value, coords.dec.value, 50/3.6e3))

    print("... Done")

def mosaic_continuum(cfg_par):
        
        os.chdir(cfg_par['general']['workdir'])

        print('# Convert continuum images to miriad format')
    
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

        #   mireg = os.path.basename(mirlist[i])
        #   mireg = string.split(mireg,'.')[0]
        #   mireg = str(i)+'_reg.mir'           

        #   miregrid = lib.miriad('regrid')
        #   miregrid.tin = '../'+baseimage_mir
        #   miregrid.in_ = mirlist[i]
        #   miregrid.bw = cfg_par[key].get('mosaic_bw',0)
        #   miregrid.out = mireg
        #   miregrid.go(rmfiles=True)
        #   mireg_list.append(mireg)
        #   puthd = lib.miriad('puthd')
        #   puthd.in_=mireg+'/pbfwhm'
        #   puthd.value=3409.5199738
        #   puthd.go()
        #   puthd = lib.miriad('puthd')

        #   puthd.type = 'double'
        #   puthd.in_=mireg+'/pntra'
        #   puthd.value=ra_deg
        #   puthd.type = 'r'
        #   puthd.go()
        #   puthd = lib.miriad('puthd')

        #   puthd.in_=mireg+'/pntdec'
        #   puthd.value=dec_deg
        #   puthd.type = 'r'
        #   puthd.go()
        #   puthd = lib.miriad('puthd')

        #   puthd.in_=mireg+'/pbinit'
        #   puthd.value=1.370031054688e9
        #   puthd.type = 'double'
        #   puthd.go()
            
        

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
    """Finds the continuum sources according to the options set in the source_finder sub-keys

    The function stores the sources as a raw txt file straight from miriad, a formatted csv file and a karma annotation file.

    Parameters
    ----------
    cfg_par : str
        Parameter file loaded with sharpener
    
    Returns
    -------
    src_list : str
        csv table of continuum sources (output of MIRIAD imsad)
    """

    print('# Find continuum sources ')

    # Getting directories and convert files if necessary
    # ++++++++++++++++++++++++++++++++++++++++++++++++++
    os.chdir(cfg_par['general']['workdir'])
    cont_im_mir = cfg_par['general']['mircontname']
    cont_im = os.path.basename(cfg_par['general']['contname'])

    # cannot use cfg_par, probably because file name would be too long for miriad
    src_imsad_out = cfg_par['general']['absdir']+'mir_src_sharpener.txt'
    # src_imsad_out = '{0:s}mir_src_sharpener.txt'.format(
    #   cfg_par['general'].get('absdir'))
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

    # Run IMSAD in Miriad to get the source
    # ++++++++++++++++++++++++++++++++++++++
    imsad = lib.miriad('imsad')
    imsad.in_ =cont_im_mir
    imsad.out = src_imsad_out
    imsad.clip = cfg_par[key]['clip']
    #imsad.region = 'boxes\('+self.abs_ex_imsad_region+'\)'
    imsad.options = cfg_par[key]['options']

    imsad.go(rmfiles=True)

    # It seems that the length of a line in the Miriad
    # output file from imsad can vary depending on the flags.
    # It is necessary to check for this and adjust the length
    check_miriad_output(src_imsad_out)
    
    
    # Modify output of imsad to save it as csv
    # ++++++++++++++++++++++++++++++++++++++++
    # changed to use ascii.read functionality
    # src_list = open(src_imsad_out,'r')
    # lines = src_list.readlines()
    # len_lines = len(lines)
    
    # ra_tmp = []
    # dec_tmp = []
    # peak_tmp = []        
    
    # for i in xrange (0,len_lines):
    #   lines[i] = lines[i].strip()
    #   tmp = lines[i].split(' ')
    #   ra_tmp.append(str(tmp[3]))
    #   dec_tmp.append(str(tmp[4]))
    #   peak_tmp.append(str(tmp[6]))
    
    # ra_tmp = np.array(ra_tmp)
    # dec_tmp = np.array(dec_tmp)
    # peak_tmp = np.array(peak_tmp)
   
    # #ID
    # ids = np.array(np.arange(1,len_lines+1,1),dtype=str)
    
    # #J2000
    # #convert ra
    # ra_vec = []
    # for i in xrange (0, len_lines):
    #   line = ra_tmp[i].split(':')
    #   last_dig = int(round(float(line[2]),0))
    #   if last_dig < 10:
    #       last_dig = '0'+str(last_dig)
    #   ra_vec.append(line[0]+line[1]+str(last_dig))
    
    # #convert dec 
    # dec_vec = []
    # dec_coord = []
    # for i in xrange (0, len_lines):
    #   line = dec_tmp[i].split(':')
    #   first_dig =  line[0][-2:]
    #   last_dig = int(round(float(line[2]),0))

    #   dec_vec.append(first_dig[1]+line[1]+str(last_dig))

    #   if line[0][-3] == '-':
    #       dec_coord.append('-'+first_dig+':'+line[1]+':'+str(last_dig)) 
    #       J2000_tmp = np.array([ a+'-'+b for a,b in zip(ra_vec,dec_vec)])
    #   if line[0][-3] == '+':
    #       dec_coord.append('+'+first_dig+':'+line[1]+':'+str(last_dig))        
    #       J2000_tmp = np.array([ a+'+'+b for a,b in zip(ra_vec,dec_vec)])

    # dec_coord = np.array(dec_coord)

    # read in data and rename columns
    src_list = ascii.read(src_imsad_out, include_names=['col3', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9', 'col10', 'col11', 'col12'])
    src_list.rename_column('col3', 'ra')
    src_list.rename_column('col4', 'dec')
    src_list.rename_column('col5', 'peak')
    src_list.rename_column('col6', 'flux_int')
    src_list.rename_column('col7', 'beam_major_decon')
    src_list.rename_column('col8', 'beam_minor_decon')
    src_list.rename_column('col9', 'beam_pang_decon')
    src_list.rename_column('col10', 'FLAG')
    src_list.rename_column('col11', 'DFLAG')
    src_list.rename_column('col12', 'FFLAG')
    n_src = np.size(src_list['ra'])

    # correct the miriad bug
    src_list['dec'] = np.array([src.replace('+0+','+') for src in src_list['dec']])

    # create two new columns
    # column 1: ID will be filled after sorting the columns
    src_id = np.zeros(n_src)

    # column 2: Source name
    j2000 = np.array(["{0:s}{1:s}".format(src_list['ra'][k].replace(':',''), src_list['dec'][k].replace(':','')) for k in range(n_src)])

    # create a table and merge it
    new_columns = Table([src_id, j2000], names=('ID', 'J2000'))
    src_list = hstack([new_columns, src_list])

    # sort the sources after source names
    src_list = src_list.group_by('J2000')

    # now assign ID
    src_list['ID'] = np.arange(n_src) + 1

    ### Find pixels of sources in continuum image
    #open continuum image
 
    # if 'CTYPE4' in prihdr:
    #   del prihdr['CTYPE4']
    # if 'CDELT4' in prihdr:
    #   del prihdr['CDELT4']    
    # if 'CRVAL4' in prihdr:
    #   del prihdr['CRVAL4']
    # if 'CUNIT4' in prihdr:
    #   del prihdr['CUNIT4']
    # if 'CRPIX4' in prihdr:            
    #   del prihdr['CRPIX4']
    # if 'CTYPE3' in prihdr:
    #   del prihdr['CTYPE3']
    # if 'CDELT3' in prihdr:
    #   del prihdr['CDELT3']
    # if 'CRVAL3' in prihdr:
    #   del prihdr['CRVAL3']
    # if 'CRPIX3' in prihdr:
    #   del prihdr['CRPIX3'] 
    # if 'CUNIT3' in prihdr:
    #   del prihdr['CUNIT3']
    # if 'NAXIS3' in prihdr:
    #   del prihdr['NAXIS3']
    # if 'NAXIS' in prihdr:
    #   del prihdr['NAXIS']
    # prihdr['NAXIS']=2
    hdulist = pyfits.open(cont_im)  # read input
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

    
    coord_list = coord.SkyCoord(src_list['ra'], src_list['dec'], unit=(u.hourangle, u.deg), frame='fk5')

    px_ra = np.zeros(n_src, dtype=int)
    px_dec = np.zeros(n_src, dtype=int)
    for k in range(n_src):
        px_ra[k], px_dec[k] = w.wcs_world2pix(coord_list[k].ra, coord_list[k].dec,1)
        px_ra[k] = int(round(px_ra[k],0))
        px_dec[k] = int(round(px_dec[k],0))

    # create a table and merge it
    new_columns = Table([px_ra, px_dec], names=('pixel_ra', 'pixel_dec'))
    src_list = hstack([src_list, new_columns])

    # print(src_list)
    # print(cfg_par['general'].get('absdir'))
    src_list.write('{0:s}mir_src_sharpener.csv'.format(cfg_par['general'].get('absdir')), format='csv', overwrite=True)

    # create a karma annotation file
    create_karma_annotation_file(coord_list,cfg_par)

    if cfg_par[key]['plot_image']:

        abs_plot.plot_continuum(cfg_par)


    # pixels_cont=np.zeros([ra_tmp.shape[0],2])
    # ra_list_tmp = np.zeros([ra_tmp.shape[0]])
    # for i in xrange (0,ra_tmp.shape[0]):
        
    #   ra_deg_tmp = conv_units.ra2deg(ra_tmp[i])
    #   dec_deg_tmp = conv_units.dec2deg(dec_coord[i])
        
    #   px,py=w.wcs_world2pix(ra_deg_tmp,dec_deg_tmp,0)
    #   pixels_cont[i,0]= str(round(px,0))
    #   pixels_cont[i,1]= str(round(py,0))        


    # #make csv file with ID, J2000, Ra, Dec, Pix_y, Pix_y, Peak[Jy] 
    # tot = np.column_stack((ids,J2000_tmp,ra_tmp,dec_coord,
    #                      pixels_cont[:,0],pixels_cont[:,1],peak_tmp))

    # write_src_csv(tot,cfg_par)  
    
    hdulist.close()

    print('# Continuum sources found #')  

    
    return src_list  






