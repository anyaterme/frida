import numpy as np
from scipy.interpolate import interp1d,interp2d
import astropy.units as u
from frida.aux_functions import interpol2newwave

doublepi = 2. * np.pi

def extrap(x, xp, yp):
    print('x=',x)
    type(x)
    y = np.interp(x, xp, yp)
    y[x < xp[0]] = yp[0] + (x[x < xp[0]] - xp[0]) * (yp[0] - yp[1]) / (xp[0] - xp[1])
    y[x > xp[-1]] = yp[-1] + (x[x > xp[-1]] - xp[-1]) * (yp[-1] - yp[-2]) / (xp[-1] - xp[-2])
    return y

def compute_r0(wave_in,seeing_ref,wave_ref=0.5 * u.micron):
    """
    Compute Fried parameter at any wavelength, given seeing at a reference wavelength
    :param in_lambda: Input wavelength
    :param seeing_ref_lambda: Seeing at the reference wavelength
    :param ref_lambda: Reference wavelength
    :return:
    """
    r0_ref = 0.98*wave_ref/seeing_ref.to('radian')*u.radian ## r_0 computed from seeing at reference lambda [in cm]
    r0 = r0_ref * (wave_in/wave_ref)**(6./5)
    return r0

def compute_seeing_lambda(wave_in,seeing_wave_ref,wave_ref= 0.5*u.micron):
    seeing = seeing_wave_ref * (wave_ref/wave_in) **(1./5)
    return seeing.to('arcsec')

def build_psf2d_2gauss(psf_2gauss,pixscale,Nx=2048,Ny=2048):
    """
    Produce a 2D image with a centred PSF, using 2048x2048 pixels with the 
    a given pixel scale. 
    It is modelled as the sum of 2 Gaussian functions: one for the core 
    (width depends linearly on the wavelength) plus one for halo (~seeing)
    :param strehl: Strehl ratio
    :param sigma_core: Width of the core Gaussian (diffraction limit core)
    :param sigma_halo: Width of the halo Gaussian (seeing)
    :return:
        
    """
    im_psf = np.array([np.zeros(Nx),np.zeros(Ny)])
    
    center = np.array([Nx/2.,Ny/2.])-0.5
    x = (np.arange(Nx)-Nx/2+0.5)*pixscale
    y = (np.arange(Ny)-Ny/2+0.5)*pixscale
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    rho = np.sqrt(xv*xv+yv*yv)
    
    amp_core = psf_2gauss["Amp_core"]
    amp_halo = psf_2gauss["Amp_halo"]
    sigma_core = psf_2gauss["FWHM_core"]/2.35
    sigma_halo = psf_2gauss["FWHM_halo"]/2.35
    
    halo2d = amp_halo * np.exp(-rho**2/2/sigma_halo**2)
    core2d = amp_core * np.exp(-rho**2/2/sigma_core**2)

    psf2d  = halo2d+core2d

    return psf2d



class GTC_AO:
    """
    This class will contain all calculations relative to GTC-AO corrections
    """

    def __init__(self,sky_conditions,guide_star, diam_telescope=10. * u.m,\
                 v0= 10 * u.m/u.s):
        ##sky_conditions={'seeing':0.9,'airmass':1.2}
        self.seeing_input = sky_conditions['seeing']
        self.wave_seeing_input = sky_conditions['wave_seeing']
        #self.airmass = sky_conditions['airmass']
        ##guide_star={'MagnitudeR':8,'Separation':0.}
        #self.mag_gs = guide_star['magnitudeR']
        #self.sep_gs = guide_star['separation']
        self.mag_gs = guide_star.Magnitude
        self.sep_gs = guide_star.Separation

        self.diam_telescope = diam_telescope
        self.v0 = v0


        # compute r0 parameter at visible range, 5000A
        self.wave_ref = 0.5 * u.micron
        self.r0_ref = compute_r0(self.wave_ref,self.seeing_input,self.wave_seeing_input)
        print ("self.r0_ref=",self.r0_ref)
        theta_isoplan_ref = 2.0 * u.arcsec # arcsec at 0.5 micron (pag 43 GTCAO System Errors Budgets 3.A)
        ## correction of isoplanatic angle due to finite aperture Keck Report 208, Expression 4.6-3, pag 4.52
        aux = diam_telescope / self.r0_ref
        self.theta_isoplan_ref_eff =  theta_isoplan_ref*(1+1/np.log10(aux))


    def compute_strehl(self,wave_in,airmass):
        #def compute_strehl(self,lamb_mic,airmass):
        # input parameter units
        #   diam_telescope in meters
        #   lambda_ref_seeing in microns

        # define required parameters


        ## determine integration time according to brightness of guide star
        mag_gs_tmin = 12.5
        c2_texp = 0.9
        if (self.mag_gs < mag_gs_tmin):
            texp_ms = 2.
        else:
            texp_ms = 2. + (self.mag_gs - mag_gs_tmin) ** 2 * c2_texp
        texp_ms *= u.ms     

        ## Parameters p0, p1 & p2 are extracted from table 1.3 RPT/OPTI/0252-R (GTCAO System Error Budgets) 3.A
        ##if (filter == 'H'):
        ##    lamb_mic = 1.65
        ##    noise_p0 = 0.2852
        ##    noise_p1 = 8.8
        ##    noise_p2 = 32.
        ##elif (filter == 'J'):
        ##    lamb_mic = 1.25
        ##    noise_p0 = 0.4216
        ##    noise_p1 = 16.5
        ##    noise_p2 = 24.
        ##elif (filter == 'K'):
        ##    lamb_mic = 2.2
        ##    noise_p0 = 0.1978
        ##    noise_p1 = 4.8
        ##    noise_p2 = 25.

        noise_tab_wave = [1.25, 1.65, 2.2] * u.micron
        noise_tab_p0 = [0.4216, 0.2852, 0.1978]
        noise_tab_p1 = [16.5, 8.8, 4.8]
        noise_tab_p2 = [24., 32. , 25.]

        noise_p0 = interpol2newwave(noise_tab_p0,noise_tab_wave,wave_in)
        noise_p1 = interpol2newwave(noise_tab_p1,noise_tab_wave,wave_in)
        noise_p2 = interpol2newwave(noise_tab_p2,noise_tab_wave,wave_in)


        #r0_cm = self.r0_vis_cm * (lamb_mic / self.lambda_vis_mic) ** (6. / 5.)
        lamb_nm = wave_in
        wvnumber_in = doublepi / wave_in ## * u.rad
        # compute r0 at visible and observing wavelengths
        r0 = compute_r0(wave_in,self.seeing_input,self.wave_seeing_input)

        ## First include static terms, which do not depend on seeing and Guide Star
        ## Static
        sigma_Calib_NCP = 25 * u.nm ## Calibration of Non-Common Path Errors
        sigma_Drift_NCP = 10 * u.nm  ## Drift of Non-Common Path Errors
        sigma_DM_hysteresis = 10 * u.nm  ## Deformable Mirror Hysteresis
        sigma_Alig_DM_WFS = 40  * u.nm ## Aligmnent of Deformable Mirror with WFS
        sigma_Calib_spot = 10  * u.nm ##
        sigma_InterMatrix = 30 * u.nm  ##
        sigma_AO_resid = 50. * u.nm  ## High order residuals from AO system ()
        sigma_hs_instrument = 50. * u.nm  ## Instrument high-spatial frequency errors
        sigma_tel_resid = 90 * u.nm  ## telescope residual optics
        sigma_tel_segvib = 60 * u.nm  ## telescope segment vibration  _
        sigma_static = np.sqrt(sigma_Calib_NCP ** 2 + sigma_Drift_NCP ** 2 + sigma_DM_hysteresis ** 2 + \
            sigma_Alig_DM_WFS ** 2 + sigma_Calib_spot ** 2 + sigma_InterMatrix ** 2 + \
            sigma_AO_resid ** 2 + sigma_hs_instrument ** 2 + sigma_tel_resid ** 2 + sigma_tel_segvib ** 2)  ##
        sigma2_static = (sigma_static * wvnumber_in) ** 2

        ####################################################
        ## fitting error (2.1.2.1 - Ref 1)
        ## It should be 119nm and 274nm when seeing/r0=[0.5"/20cm] and [1.5"/6.9cm], 
        ## respectively (Data taken from Table 3 - Ref 1)
        ##  aproximation log(RMS_nm) = 3.0956 -0.784 * log10(r0@0.5mic)
        r0_0p5mic_cm = self.r0_ref.to("cm").value
        rms_fitting_nm = 10. ** (3.0956 - 0.784 * np.log10(r0_0p5mic_cm)) * u.nm
        sigma2_fitting = (rms_fitting_nm * wvnumber_in) ** 2
        ## includes the variation with airmass, which affects sources of error depending on seeing, or r0
        sigma2_fitting *= airmass

        ####################################################
        ## bandwidth errors (taken from Table 5 - Ref 1)
        ##  value for turbulent speed (\nu) = 10 m/s
        ##   r0 [cm]     |  20    |  15.6  | 6.7
        ##   Sigma [nm]  |   54.  |  67.   | 135
        r0_bandwidth_tab = [20., 15.6, 6.7]* u.cm
        v0_bandwidth_tab = [10.,15.] * u.m / u.s
        sigma_bandwidth_tab = [[54., 67., 135.],[76,94,189.]] * u.nm

        r0_interp = r0_0p5mic_cm
        v0_interp = self.v0.to('m/s')
        _zunit = sigma_bandwidth_tab.unit
        _xx = r0_bandwidth_tab.to("cm").value
        _yy = v0_bandwidth_tab.to("m/s").value
        _zz = sigma_bandwidth_tab.to("nm").value
        f = interp2d(_xx, _yy, _zz, kind='linear')
        sigma_bandwidth = f(r0_interp,v0_interp) 
        sigma_bandwidth *= _zunit
        sigma2_bandwidth = (sigma_bandwidth * wvnumber_in) ** 2
        ## includes the variation with airmass, which affects sources of error depending on seeing, or r0
        ## Section 2.4, pag 38 ref 1
        sigma2_bandwidth *= airmass


        ####################################################
        ## Photon noise error
        effic_wfs = 0.22  ## (pag 34, Ref 1)
        rad_subaper = self.diam_telescope / 20. / 2.  ## diameter of primary divided into 20 subapertures
        area_subaper = (rad_subaper) ** 2. * np.pi  ## area of subaperture in m^2
        nphotR0 = 3.0216e7 * u.photon / u.m**2 / u.ms ## photons/msec/m2 mag_R=0.
        nphot = nphotR0 * 10 ** (-0.4 * self.mag_gs) * texp_ms * area_subaper * effic_wfs
        sigma2_photnoise = noise_p1 / nphot
        sigma2_readnoise = noise_p2 / nphot / nphot
        # adjust units 
        sigma2_photnoise *= u.photon
        sigma2_readnoise *= u.photon**2


        ####################################################
        ## Scintillation
        h_scint = 5 * u.km
        sigma2_scint = 0.288 * (np.sqrt(wave_in * h_scint) / r0) ** (5. / 3.)


        ####################################################
        ## low-order tip-tilt or bandwidth low order (it uses a different expression)
        sigma_tt_bw_rad = 11.3 * 1.e-9 ## due to bandwidth (specified in 1E-9 rad)
        sigma_tt_wndsh_rad = 31.8 * 1.e-9  ## due to wind shake 0.1 arcsecond
        sigma_tt_rad = np.sqrt(sigma_tt_bw_rad**2+sigma_tt_wndsh_rad**2)
        lamb_over_D =  wave_in / self.diam_telescope
        ## to compute the ratio (sigma_nrad*1e-9) / (lamb_nm*1.e-9/diam_telescope), the factor 1.e-9 cancelled out
        strehl_tt = 1. / (1 + np.pi * np.pi / 2. * (sigma_tt_rad / lamb_over_D) ** 2)


        ## error due to separation of guide star relative to the pointing
        ##sigma_anisop = (dist_gs/theta_isoplan_ref)**(5./3)*(lamb_ref/lamb_mic)**2
        theta_isoplan_eff = self.theta_isoplan_ref_eff * (wave_in/self.wave_ref) ** 6/5
        sigma2_anisop = (self.sep_gs / theta_isoplan_eff) ** (5. / 3.) ## * (self.lambda_vis_mic / lamb_mic) ** 2

        ## adding all contributions computed before
        ###sigma2_total = sigma2_static + sigma2_scint + sigma2_bandwidth + sigma2_fitting + \
           ###        sigma2_photnoise  + sigma2_anisop  + sigma2_readnoise
        print("sig_static=",np.sqrt(sigma2_static+0.))
        print("sig_scint=",np.sqrt(sigma2_scint+0.))
        print("sig_bandwidth=",np.sqrt(sigma2_bandwidth+0.))
        print("sig_fitting=",np.sqrt(sigma2_fitting+0.))
        print("sig_photnoise=",np.sqrt(sigma2_photnoise+0.))
        print("sig_readnoise=",np.sqrt(sigma2_readnoise+0.))
        print("sig_anisop=",np.sqrt(sigma2_anisop+0.))
        sigma2_total = sigma2_static + sigma2_scint + sigma2_bandwidth + sigma2_fitting + \
                   sigma2_photnoise  + sigma2_readnoise + sigma2_anisop 
        print("sig2_total=",np.sqrt(sigma2_total+0.))
        strehl = np.exp(-sigma2_total) * strehl_tt
        print("Strehl=",strehl)

        ## as output produces Strehl ratio, RMS wavefront error components in nanometers
        return {'StrehlR': strehl, 'EffIsoplanAngle':theta_isoplan_eff,'sigma2_scint': sigma2_scint,  \
            'sigma2_stat': sigma2_static,'sigma2_bandwidth': sigma2_bandwidth,  \
            'strehl_tt': strehl_tt, 'sigma2_fitting': sigma2_fitting, \
            'sigma2_photnoise': sigma2_photnoise,'sigma2_readnoise': sigma2_readnoise, \
            'sigma2_anisop': sigma2_anisop,'nphot':nphot,'Texp_ms':texp_ms,\
            'wvnumber':wvnumber_in}


    def compute_ee(self, psf, pixscale, fcore=1.5, spaxel=False,\
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

    def compute_ee_box(self, psf,fcore=1.5,slit_width=2):
        """
        :param strehl:
        :param seeing:
        :param filter:
        :param teldiam:
        :param fcore:
        :return:
        """
        ### HAY QUE ESCRIBIR LAS FORMULAS CORRECTAMENTE PARA UNA APERTURA TIPO RECTANGULAR
        raper_ref = fcore * psf["FWHM_core"]

        # raper_core=raper_ref/sigma_core

        raper2_ref = raper_ref * raper_ref
        sigma2_halo = psf["sigma_halo"] * psf["sigma_halo"]
        sigma2_core = psf["sigma_core"] * psf["sigma_core"]

        #psf = compute_psf(strehl, sigma_core, sigma_halo)

        A_core = psf['Amp_core']
        A_halo = psf['Amp_halo']
        fcore = psf['rat_core2halo']

        ee = (A_halo * sigma2_halo * (1 - np.exp(-raper2_ref / 2. / sigma2_halo)) + \
              A_core * sigma2_core * (1 - np.exp(-raper2_ref / 2. / sigma2_core))) / \
             (A_halo * sigma2_halo + A_core * sigma2_core)

        ee1 = 1. - 1. / (sigma2_core * fcore + sigma2_halo) * (sigma2_halo * np.exp(-raper2_ref / 2. / sigma2_halo) + \
                                                               fcore * sigma2_core * np.exp(-raper2_ref / 2. / sigma2_core))

        return {'EE': ee,'EE1': ee1, 'ApertRad':raper_ref}

    def compute_psf(self,wave,strehl):
        """
        The PSF is normalized to have area equals unity.
        :param strehl: Strehl ratio
        :param sigma_core: Width of the core Gaussian (diffraction limit core)
        :param sigma_halo: Width of the halo Gaussian (seeing)
        :return:
            
        """
        
        fwhm_core = 0.0242 * u.arcsec * (wave/u.micron) * (10.4 * u.m/ self.diam_telescope)
        ## seeing is assumed to be measured at 5000AA = 0.5 microns
        fwhm_halo = compute_seeing_lambda(wave,self.seeing_input,self.wave_seeing_input)

        # homogeneize units to arcsec 
        fwhm_core = fwhm_core.to('arcsec')
        fwhm_halo = fwhm_halo.to('arcsec')
  
        sigma2_halo = (fwhm_halo / 2.35)**2
        sigma2_core = (fwhm_core / 2.35)**2

        amp_core_norm = 1. / doublepi / (sigma2_halo - sigma2_core) * \
            (strehl*sigma2_halo/sigma2_core -1.)

        amp_halo_norm = 1. / doublepi / sigma2_halo * (1. - (strehl*sigma2_halo-sigma2_core)/\
                                         (sigma2_halo - sigma2_core))

        return {"Model_PSF":'2-Gaussians',\
            "Amp_core": amp_core_norm, "Amp_halo": amp_halo_norm, \
            "sigma2_core":sigma2_core,"sigma2_halo":sigma2_halo,\
            "FWHM_core":fwhm_core.to('arcsec'),"FWHM_halo":fwhm_halo.to('arcsec')}

    def set_aperture_circular(self,pixscale,radius,spaxel=False):

        area_pixel = pixscale*pixscale
        if (spaxel): area_pixel *=  2.
        area_apert=np.pi*radius**2 ## in arcsec
        npix_apert = area_apert / area_pixel

        return {'radius':radius,'area_aperture':area_aperture,'Npix':npix_apert,'area_pixel':area_pixel}
