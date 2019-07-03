import sys
import string
import os
import math
import numpy as np
#import Cosmo as c
from kk import *
from astropy import units
from astropy.io import fits


class hi:
    '''Tools to analyze the HI spectral line

    - hiline 
            frequency of the HI given a redshift
    - tau_abs
            converts absorbed flux in optical depth
    - nhi_abs
            determines the column density from optical depth
            assuming temperature TSPIN = 100
    - beam_area
            determines the area of the beam of the observations
    - mhi_abs
            mass of the absorbed HI
    - nhi_em
            column density from emission line (and beam of observation)
    - mhi_em
            mass of HI from emission line
    - mhi_flux
            mass of HI from integrated flux and redshift of source
    '''

    def __init__(self):

        self.hi = kk.HI
        self.c = kk.C
        self.knhi = kk.KnhiABS
        self.nhiem = kk.KnhiEM
        self.T = kk.TSPIN
        self.mhi = kk.MHI
        self.msun = kk.MSUN
        self.chi = kk.CHI
        self.mp = kk.MP

    def hiline(self, z):
        '''
        Estimates the expected velocity of the HI line given the redshift
        of the source

        Parameters:
            z: redshift (float)

        Returns:
            hi.velhi: velocity in km/s
            hi.freq: frequency of the HI line in MHz
        '''

        freq = self.hi/(1+z)/1e06  # MHz
        velocity = self.C*((self.hi-freq)/freq)/1e5  # km/s

        print('HI expected frequency = '+str(round(freq, 3))+' MHz')
        print('HI systemic velocity = '+str(round(velocity, 3))+' km/s')

        return freq, velocity

    def optical_depth(self, scont, sabs):
        '''
        Estimates the optical depth of an absorption line

        Parameters:
            scont: continuum flux (float)
            sabs: flux of the absorbed component (negative float)
        
        Returns:
            hi.tau: optical depth of the line	
        '''

        tau = np.log(1.-(-sabs/scont))

        if tau.size == 1:
            print('Optical depth = '+str(round(tau, 3)))

        return tau

    def nhi_abs(tau, dv):
        '''Estimates the column density of the absorption line

        Parameters:
            tau: optical depth of the line (float)
            dv: width of the line in km/s
        
        Returns:
            hi.nhi_abs: column density of the absorption line in cm-2
        '''

        nhiabs = kk.knhi*kk.T*tau*dv

        print('N(HI) = '+str(round(nhiabs, 3))+' cm-2')

        return nhiabs

    def beam_area(bx, by, z):
        '''
        Estimates the area of the beam of the observations

        Parameters:
            bx: x axis of the beam in arcsec (float)
            by: y axis of the beam in arcsec (float)

        Returns:
            hi.beamarea: area of the beam in cm2	
        '''

        bxcm = c.ang2lin(bx, z)*1e6
        bycm = c.ang2lin(by, z)*1e6
        beamarea = (bxcm*bycm)*(PC**2)

        print('Beam Area = '+str(round(beamarea, 3))+' cm2')

        return beamarea

    def mhi_abs(nhi_abs, area):
        '''Estimates the mass of the absorbed HI

        Parameters:
            nhi_abs: column density of the absorption line in cm-2
            area: area over which the column density is integrated in cm (output of hi.beam_area)
        
        Returns:
            hi.mhi_abs: hi mass inferred by the absorption line in Msun
        '''

        mhiabs = area*self.mp*nhi_abs/self.msun

        print('M(HI) = '+str(round(mhiabs, 3))+' Msun')

        return mhiabs

    def nhi_em(self, s, bx, by, dv):
        '''Estimates the column density of the absorption line

        Parameters:
            s: flux of the emission line in mJy
            bx: x axis of the beam in arcsec (float)
            by: y axis of the beam in arcsec (float)
            dv: width of the line in km/s
        
        Returns:
            hi.nhi_abs: column density of the absorption line in cm-2
        '''

        bxa = bx/60
        bya = by/60  # convert beam in arcmin
        nhiem = self.nhiem*dv*s/(bxa*bya)

        print('N(HI) = '+str(round(nhi_em, 3))+' cm-2')

        return nhi_em

    def mhi_em(self, nhi, area):
        '''
        Estimates the column density of the absorption line

        Parameters:
            nhi: column density of the emission line
            area:area over which the column density is integrated in cm	(output of hi.beam_area)
        
        Returns:
            hi.mhiem: mass inferred from the HI absoprtion line in Msun	
        '''

        mhiem = nhi*self.mhi*area/MSUN

        print('M(HI) = '+str(round(mhiem, 3))+' Msun')

        return mhiem

    def mhi_flux(self, z, s, bx, by, pix):
        '''
        Estimates the HI mass from the observed flux, the distance of the
        source, and the resolution of the observations
        (beam and pixel size in arcsec)

        Parameters:
            z: redshift of the source (float)
            bx: x axis of the beam in arcsec (float)
            by: y axis of the beam in arcsec (float)
            pix: pixel size in arcsec (float)
        
        Returns:
            hi.mhflux: mass inferred from the flux of the HI emission line	
        '''

        dl = c.lum_dist(z)/3.085678e24
        beamcorr = pix**2/(bx*by)
        bla = s*beamcorr
        mhi = self.chi*(dl**2)*bla

        return mhi
