__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"


import os, sys, string
import numpy as np
from prettytable import PrettyTable
from astropy.io import ascii



class kk:
	'''
	Useful constants called by other radiobs classes as kk.constant
	'''

	#converstions
	RAD2DEG=180./np.pi
	PC=3.08567758E18    #cm
	JANSKY=1e-23        #erg/scm2Hz
	
	#constants
	HI=1.42040575177e+09 #Hz
	TSPIN=100            #K
	C=2.99792458E10     #cm/s
	G=6.6742E-08        #cm3kg-1s-1 
	MSUN=1.98855e33      #g
	MHI=1.6749E-24       #g
	CHI=2.36E5
	MP=1.67492728E-24   #g
	SIGMAT=6.66524E-25  #cm2
	
	#convert to column density of HI
	KnhiABS = 1.8216E18
	KnhiEM = 3.1E17
	
	#magnitudes & factors for Tully-Fisher
	M26 = -23.33
	STULLY = -9.64

	#cosmological constants
	h0 = 70
	omega_l = 0.7 
	omega_m = 0.3

	#HI
	MILKY_LEFT = -30.
	MILKY_RIGHT = +30.


	def tablefy(self, tablename, columns, columnames):
		'''
		Print pretty table out of array of columnames and columns
		INPUT
			tablename : name of outputable
			column: array with data of columns, shape (len(columnames),data)
			columnames: list of names of columns
		OUTPUT
			x: table in txt format in tablename file
		RETURN
			x: table in PrettyTable format

		'''	
		x = PrettyTable()


		for i in xrange(0,len(columnames)):

			x.add_column(columnames[i],columns[i,:])

		basename = string.split(tablename,'.')
		datatmp = x.get_html_string()
		htmltable = basename[0]+'.html'
		with open(htmltable, 'wb') as f:
			f.write(datatmp)

		datatmp = ascii.read(htmltable, format='html')
		ascii.write(datatmp,tablename,format='csv',names = datatmp.dtype.names)
		os.remove(htmltable)
		print(x)

		return x