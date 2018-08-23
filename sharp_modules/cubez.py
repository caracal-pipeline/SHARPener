__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"


from astropy.io import fits
import numpy as np

def zaxis(cubename):

	cubefile = fits.open(cubename)  # read input
	hdr = cubefile[0].header
			
	freq = (np.linspace(1, hdr['NAXIS3'], hdr['NAXIS3']) - hdr['CRPIX3']) * hdr['CDELT3'] + hdr['CRVAL3']

	return freq 