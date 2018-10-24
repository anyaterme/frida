#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 18:34:11 2018

@author: JoseAcosta
"""

import numpy as np
from scipy.special import jv
import astropy.units as u

doublepi = 2. * np.pi


def compute_ee_2Gauss(psf, pixscale, fcore=1.5, spaxel=False,\
                   source_type='Point source'):
        """
        :param strehl:
        :param seeing:
        :param filter:i
        :param teldiam:
        :param fcore:
        :return:
        """
        if (source_type == "Point source"):
           raper_ref = fcore * psf["FWHM_core"]
        else:
           raper_ref = fcore * pixscale
            
        # raper_core=raper_ref/sigma_core

        raper2_ref = raper_ref**2
        sigma2_halo = psf["sigma2_halo"]
        sigma2_core = psf["sigma2_core"]
        print("=========================")
        print("raper2_ref=",raper2_ref) 
        print("sigma2_halo=",sigma2_halo) 
        print("sigma2_core=",sigma2_core) 
        print("=========================") 

        #psf = compute_psf(strehl, sigma_core, sigma_halo)

        if (source_type == "Point source"):
            A_core = psf['Amp_core']
            A_halo = psf['Amp_halo']
            ee_core = doublepi * A_core * sigma2_core * \
                 (1 - np.exp(-raper2_ref/2./sigma2_core))        
            ee_halo = doublepi * A_halo * sigma2_halo * \
                 (1 - np.exp(-raper2_ref/2./sigma2_halo))             
            ee = ee_core + ee_halo 
        else: 
            ee_core = 1.              
            ee_halo = 0.              
            ee = ee_core + ee_halo 
            ## extended source

        area_per_pixel = pixscale*pixscale
        ## check if spaxel is set, then assume a spaxel as 2x1 pixel, according to FRIDA specs.
        if (spaxel): area_per_pixel *= 2.

        area_apert=np.pi*raper_ref**2 ## in arcsec
        npix_apert = area_apert / area_per_pixel

        return {'EE': ee,'EE-core': ee_core,'EE-halo': ee_halo, 'Radius':raper_ref, \
                'Npix':npix_apert, 'Area':area_apert,\
                'Area_pixel':area_per_pixel}

def compute_ee_AiryGauss(psf, pixscale, fcore=1.5, spaxel=False,\
                   source_type='Point source'):
        """
        :param strehl:
        :param seeing:
        :param filter:i
        :param teldiam:
        :param fcore:
        :return:
        """
        if (source_type == "Point source"):
           raper_ref = fcore * psf["FWHM_core"]
        else:
           raper_ref = fcore * pixscale
            
        # raper_core=raper_ref/sigma_core

        raper2_ref = raper_ref**2
        sigma2_halo = psf["sigma2_halo"]
        print("=========================")
        print("raper2_ref=",raper2_ref) 
        print("sigma2_halo=",sigma2_halo) 
        print("=========================") 

        #psf = compute_psf(strehl, sigma_core, sigma_halo)

        if (source_type == "Point source"):
            A_core = psf['Amp_core']
            A_halo = psf['Amp_halo']
            xscale = 2 * np.pi * psf["Diameter"].to(u.meter) / psf["Wave"].to(u.meter)
            x = (xscale * raper_ref(u.radian)).value
            ee_core = A_core * \
                 (1 - jv(0,x)**2 - jv(1,x)**2)        
            ee_halo = doublepi * A_halo * sigma2_halo * \
                 (1 - np.exp(-raper2_ref/2./sigma2_halo))             
            ee = ee_core + ee_halo 
        else: 
            ee_core = 1.              
            ee_halo = 0.              
            ee = ee_core + ee_halo 
            ## extended source

        area_per_pixel = pixscale*pixscale
        ## check if spaxel is set, then assume a spaxel as 2x1 pixel, according to FRIDA specs.
        if (spaxel): area_per_pixel *= 2.

        area_apert=np.pi*raper_ref**2 ## in arcsec
        npix_apert = area_apert / area_per_pixel

        return {'EE': ee,'EE-core': ee_core,'EE-halo': ee_halo, 'Radius':raper_ref, \
                'Npix':npix_apert, 'Area':area_apert,\
                'Area_pixel':area_per_pixel}


def buildim_psf_2gauss(psf_2gauss,pixscale,Nx=250,Ny=250):
    """
    Produce a 2D image with a centred PSF, using 250x250 pixels with the 
    a given pixel scale. 
    It is modelled as the sum of 2 Gaussian functions: one for the core 
    (width depends linearly on the wavelength) plus one for halo (~seeing)
    :param strehl: Strehl ratio
    :param sigma_core: Width of the core Gaussian (diffraction limit core)
    :param sigma_halo: Width of the halo Gaussian (seeing)
    :return:
        
    """
    x = (np.arange(Nx)-Nx/2+0.5)*pixscale
    y = (np.arange(Ny)-Ny/2+0.5)*pixscale
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    rho = np.sqrt(xv*xv+yv*yv)
    
    amp_core = psf_2gauss["Amp_core"]
    amp_halo = psf_2gauss["Amp_halo"]
    sigma_core = psf_2gauss["FWHM_core"]/2.35
    sigma_halo = psf_2gauss["FWHM_halo"]/2.35
    
    print("Amp halo",amp_halo)
    print("Amp core",amp_core)
    halo2d = amp_halo * np.exp(-rho**2/2/sigma_halo**2)
    core2d = amp_core * np.exp(-rho**2/2/sigma_core**2)

    psf2d  = halo2d+core2d

    return psf2d * pixscale * pixscale,x,y


def buildcube_psf_2gauss(psf_wave,pixscale,Nx=250,Ny=250):
    """
    Produce a 2D image with a centred PSF, using 250x250 pixels with the 
    a given pixel scale. 
    It is modelled as the sum of 2 Gaussian functions: one for the core 
    (width depends linearly on the wavelength) plus one for halo (~seeing)
    :param strehl: Strehl ratio
    :param sigma_core: Width of the core Gaussian (diffraction limit core)
    :param sigma_halo: Width of the halo Gaussian (seeing)
    :return:
        
    """
    Nwave = len(psf_wave['Amp_core'])
    
    psf_cube = np.ndarray(shape=(Nwave,Nx,Ny),dtype='float')
    
    x = (np.arange(Nx)-Nx/2+0.5)*pixscale
    y = (np.arange(Ny)-Ny/2+0.5)*pixscale
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    rho = np.sqrt(xv*xv+yv*yv)
    
    for i in range(Nwave):
        amp_core = psf_wave["Amp_core"][i]
        amp_halo = psf_wave["Amp_halo"][i]
        sigma_core = psf_wave["FWHM_core"][i]/2.35
        sigma_halo = psf_wave["FWHM_halo"][i]/2.35
    
        halo2d = amp_halo * np.exp(-rho**2/2/sigma_halo**2)
        core2d = amp_core * np.exp(-rho**2/2/sigma_core**2)

        psf_cube[i,:,:]  = halo2d+core2d

    print('psf_cube[wave]',psf_cube[900:1000,int(Nx/2),int(Ny/2)]*pixscale**2)
    print('plano suma',psf_cube[950,:,:].sum()*pixscale**2)
    print('psf_cube[row]',psf_cube[950,:,int(Ny/2)]*pixscale**2)

    return psf_cube * pixscale * pixscale,x,y


def buildim_psf_AiryGauss(psf,pixscale,Nx=250,Ny=250):
    """
    Produce a 2D image with a centred PSF, using Nx x Ny pixels with the 
    a given pixel scale. 
    It is modelled as the sum of 2 Gaussian functions: one for the core 
    (width depends linearly on the wavelength) plus one for halo (~seeing)
    :param psf: Definition of psf with par
    :param pixscale: Pixel scale of the output image 
    :
    :return:
        
    """
    
    x = (np.arange(Nx)-Nx/2+0.5)*pixscale
    y = (np.arange(Ny)-Ny/2+0.5)*pixscale
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    rho = np.sqrt(xv*xv+yv*yv)
    
    amp_core = psf["Amp_core"]
    amp_halo = psf["Amp_halo"]
    sigma_halo = psf["FWHM_halo"]/2.35
    
    halo2d = amp_halo * np.exp(-rho**2/2/sigma_halo**2)
    #print("Diameter",psf["Diameter"])
    #print("Wave",psf["Wave"])
    xscale = 2 * np.pi * psf["Diameter"].to(u.meter) / psf["Wave"].to(u.meter)
    #print("xscale ",xscale)
    #print("rho ",rho)
    x = (xscale * rho.to(u.radian)).value
    #print("check units of x ")
    #print("x ",x)
    core2d = amp_core * (2*jv(1,x)/x)**2

    psf2d  = halo2d+core2d

    return psf2d * pixscale * pixscale,x,y


def buildcube_psf_AiryGauss(psf_wave,pixscale,Nx=250,Ny=250):
    """
    Produce a 2D image with a centred PSF, using Nx x Ny pixels with the 
    a given pixel scale. 
    It is modelled as the sum of 2 Gaussian functions: one for the core 
    (width depends linearly on the wavelength) plus one for halo (~seeing)
    :param psf: Definition of psf with par
    :param pixscale: Pixel scale of the output image 
    :
    :return:
        
    """
    Nwave = len(psf_wave['Wave'])
    
    psf_cube = np.ndarray(shape=(Nwave,Nx,Ny),dtype='float')
    
    x = (np.arange(Nx)-Nx/2+0.5)*pixscale
    y = (np.arange(Ny)-Ny/2+0.5)*pixscale
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    rho = np.sqrt(xv*xv+yv*yv)
    
    for i in range(Nwave):
        amp_core = psf_wave["Amp_core"][i]
        amp_halo = psf_wave["Amp_halo"][i]
        sigma_halo = psf_wave["FWHM_halo"][i]/2.35
    
        halo2d = amp_halo * np.exp(-rho**2/2/sigma_halo**2)
        #print("Diameter",psf_wave["Diameter"])
        #print("Wave",psf_wave["Wave"][i])
        xscale = 2 * np.pi * psf_wave["Diameter"].to(u.meter) / \
             psf_wave["Wave"][i].to(u.meter)
        #print("xscale ",xscale)
        #print("rho ",rho)
        x = (xscale * rho.to(u.radian)).value
        #print("check units of x ")
        #print("x ",x)
        core2d = amp_core * (2*jv(1,x)/x)**2

        psf_cube[i,:,:]  = halo2d+core2d
        
    print('psf_cube[wave]',psf_cube[900:1000,int(Nx/2),int(Ny/2)]*pixscale**2)
    print('plano suma',psf_cube[950,:,:].sum()*pixscale**2)
    print('psf_cube[row]',psf_cube[950,:,int(Ny/2)]*pixscale**2)

    return psf_cube * pixscale * pixscale,x,y

