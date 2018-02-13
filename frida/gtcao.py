import numpy as np
from math import  pi
from scipy.interpolate import interp1d
import astropy.units as u


def extrap(x, xp, yp):
    print('x=',x)
    type(x)
    y = np.interp(x, xp, yp)
    y[x < xp[0]] = yp[0] + (x[x < xp[0]] - xp[0]) * (yp[0] - yp[1]) / (xp[0] - xp[1])
    y[x > xp[-1]] = yp[-1] + (x[x > xp[-1]] - xp[-1]) * (yp[-1] - yp[-2]) / (xp[-1] - xp[-2])
    return y


def compute_r0(in_lambda,seeing_ref_lambda,ref_lambda=0.5):
    """
    Compute Fried parameter at any wavelength, given seeing at a reference wavelength
    :param in_lambda: Input wavelength
    :param seeing_ref_lambda: Seeing at the reference wavelength
    :param ref_lambda: Reference wavelength
    :return:
    """
    r0_ref_lambda = 20.637*ref_lambda/seeing_ref_lambda ## r_0 computed from seeing at reference lambda [in cm]
    lambda_ref_vis = 0.5 ## Reference lambda set at 5000A=0.5micras
    seeing_ref_vis = 0.5  ## Seeing arcseconds at reference lambda 0.5 microns
    r0_ref_vis = 20.  ## r0 in cm equivalent to seeing_ref_vis=0.5 arcsec at 0.5microns
    r0_cm = r0_ref_lambda * (in_lambda/ref_lambda)**(6./5)
    return r0_cm.value * u.cm

def compute_seeing_lambda(in_lambda,seeing_ref_lambda,ref_lambda=0.5):
    seeing = seeing_ref_lambda * (ref_lambda/in_lambda) **(1./5)
    return seeing

class GTC_AO:
    """
    This class will contain all calculations relative to GTC-AO corrections
    """

    def __init__(self,sky_conditions,guide_star, diam_telescope=10.):
        ##sky_conditions={'seeing':0.9,'airmass':1.2}
        self.seeing_input = sky_conditions['seeing']
        self.lambda_seeing_input = sky_conditions['lambda_seeing']
        #self.airmass = sky_conditions['airmass']
        ##guide_star={'MagnitudeR':8,'Separation':0.}
        #self.mag_gs = guide_star['magnitudeR']
        #self.sep_gs = guide_star['separation']
        self.mag_gs = guide_star.Magnitude
        self.sep_gs = guide_star.Separation

        self.diam_telescope = diam_telescope


        # compute r0 parameter at visible range, 5000A
        self.lambda_vis_mic = 0.5 * u.micron
        self.r0_vis_cm = compute_r0(self.lambda_vis_mic,self.seeing_input,self.lambda_seeing_input)
        print ("self.r0_vis_cm=",self.r0_vis_cm)
        theta_isoplan_vis = 2.0  # arcsec at 0.5 micron (pag 43 GTCAO System Errors Budgets 3.A)
        ## correction of isoplanatic angle due to finite aperture Keck Report 208, Expression 4.6-3, pag 4.52
        aux = (diam_telescope / self.r0_vis_cm).si.value
        self.theta_isoplan_vis_eff =  theta_isoplan_vis*(1+1/np.log10(aux))


    def compute_strehl(self,lamb_mic,airmass):
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
        texp_ms = texp_ms * u.ms     

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

        noise_tab_lamb = [1.25, 1.65, 2.2]  #) * u.micron
        noise_tab_p0 = [0.4216, 0.2852, 0.1978]
        noise_tab_p1 = [16.5, 8.8, 4.8]
        noise_tab_p2 = [24., 32. , 25.]

        lamb_interp = lamb_mic.to("micron").value
        #noise_p0 = extrap(lamb_mic, noise_tab_lamb, noise_tab_p0)[0]
        noise_p0 = np.interp(lamb_interp,noise_tab_lamb,noise_tab_p0)
        noise_p1 = np.interp(lamb_interp, noise_tab_lamb, noise_tab_p1)
        noise_p2 = np.interp(lamb_interp, noise_tab_lamb, noise_tab_p2)


        #r0_cm = self.r0_vis_cm * (lamb_mic / self.lambda_vis_mic) ** (6. / 5.)
        lamb_nm = lamb_mic
        # compute r0 at visible and observing wavelengths
        r0_cm = compute_r0(lamb_mic,self.seeing_input,self.lambda_seeing_input)

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
        sigma_static_nm = np.sqrt(sigma_Calib_NCP ** 2 + sigma_Drift_NCP ** 2 + sigma_DM_hysteresis ** 2 + \
            sigma_Alig_DM_WFS ** 2 + sigma_Calib_spot ** 2 + sigma_InterMatrix ** 2 + \
            sigma_AO_resid ** 2 + sigma_hs_instrument ** 2 + sigma_tel_resid ** 2 + sigma_tel_segvib ** 2)  ##
        sigma2_static = (sigma_static_nm * 2. * pi / lamb_mic) ** 2


        ####################################################
        ## fitting error (2.1.2.1 - Ref 1)
        ## It should be 119nm and 274nm when seeing/r0=[0.5"/20cm] and [1.5"/6.9cm], respectively (Data taken from Table 3 - Ref 1)
        ##  aproximation log(RMS_nm) = 3.0956 -0.784 * log10(r0@0.5mic)
        r0_0p5mic_cm = self.r0_vis_cm.to("cm").value
        rms_fitting_nm = 10. ** (3.0956 - 0.784 * np.log10(r0_0p5mic_cm)) * u.nm
        sigma2_fitting = (rms_fitting_nm / lamb_nm * 2 * pi) ** 2
        ## includes the variation with airmass, which affects sources of error depending on seeing, or r0
        sigma2_fitting = sigma2_fitting * airmass

        ####################################################
        ## bandwidth errors (taken from Table 5 - Ref 1)
        ##  value for turbulent speed (\nu) = 10 m/s
        ##   r0 [cm]     |  20    |  15.6  | 6.7
        ##   Sigma [nm]  |   54.  |  67.   | 135
        r0_bandwidth_tab = [20., 15.6, 6.7]
        sigma_bandwidth_tab = [54., 67., 135.]

        r0_interp = r0_0p5mic_cm
        sigma_bandwidth_nm = 10. ** np.interp(r0_interp, r0_bandwidth_tab, sigma_bandwidth_tab) * u.nm
        sigma2_bandwidth = (sigma_bandwidth_nm * 2. * pi / lamb_nm) ** 2
        ## includes the variation with airmass, which affects sources of error depending on seeing, or r0
        ## Section 2.4, pag 38 ref 1
        sigma2_bandwidth = sigma2_bandwidth * airmass


        ####################################################
        ## Photon noise error
        effic_wfs = 0.22  ## (pag 34, Ref 1)
        rad_subaper = self.diam_telescope / 20. / 2.  ## diameter of primary divided into 20 subapertures
        area_subaper = (rad_subaper) ** 2. * pi  ## area of subaperture in m^2
        nphotR0 = 3.0216e7 / u.m**2 / u.ms ## photons/msec/m2 mag_R=0.
        nphot = nphotR0 * 10 ** (-0.4 * self.mag_gs) * texp_ms * area_subaper * effic_wfs
        sigma2_photnoise = noise_p1 / nphot
        sigma2_readnoise = noise_p2 / nphot / nphot


        ####################################################
        ## Scintillation
        h_scint = 5 * u.km
        sigma2_scint = 0.288 * (np.sqrt(lamb_mic * h_scint) / r0_cm) ** (5. / 3.)


        ####################################################
        ## low-order tip-tilt or bandwidth low order (it uses a different expression)
        sigma_tt_bw_nrad = 11.3  ## due to bandwidth (specified in 1E-9 rad)
        sigma_tt_wndsh_nrad = 31.8  ## due to wind shake 0.1 arcsecond
        sigma_tt_nrad = sigma_tt_wndsh_nrad
        sigma_tt_rad = sigma_tt_nrad * 1.e-9
        lamb_over_D = lamb_nm / self.diam_telescope
        ## to compute the ratio (sigma_nrad*1e-9) / (lamb_nm*1.e-9/diam_telescope), the factor 1.e-9 cancelled out
        strehl_tt = 1. / (1 + pi * pi / 2. * (sigma_tt_nrad / lamb_over_D) ** 2)


        ## error due to separation of guide star relative to the pointing
        ##sigma_anisop = (dist_gs/theta_isoplan_ref)**(5./3)*(lamb_ref/lamb_mic)**2
        theta_isoplan_eff = self.theta_isoplan_vis_eff * (lamb_mic/self.lambda_vis_mic) ** 6/5
        sigma2_anisop = (self.sep_gs / theta_isoplan_eff) ** (5. / 3.) ## * (self.lambda_vis_mic / lamb_mic) ** 2

        ## adding all contributions computed before
        ###sigma2_total = sigma2_static + sigma2_scint + sigma2_bandwidth + sigma2_fitting + \
           ###        sigma2_photnoise  + sigma2_anisop  + sigma2_readnoise
        print("sig2_static=",sigma2_static)
        print("sig2_scint=",sigma2_scint)
        print("sig2_bandwidth=",sigma2_bandwidth)
        print("sig2_fitting=",sigma2_fitting)
        print("sig2_photnoise=",sigma2_photnoise)
        print("sig2_readnoise=",sigma2_readnoise)
        print("sig2_anisop=",sigma2_anisop)
        sigma2_total = sigma2_static + sigma2_anisop 
        strehl = np.exp(-sigma2_total) * strehl_tt
        lamb_nm_2pi = lamb_nm / 2 / pi

        ## as output produces Strehl ratio, RMS wavefront error components in nanometers
        return {'StrehlR': strehl, 'EffIsoplanAngle':theta_isoplan_eff,'rms_scint': np.sqrt(sigma2_scint) * lamb_nm_2pi,  \
            'rms_stat': np.sqrt(sigma2_static) * lamb_nm_2pi,'rms_bandwidth': np.sqrt(sigma2_bandwidth) * lamb_nm_2pi,  \
            'strehl_tt': strehl_tt, 'rms_fitting': np.sqrt(sigma2_fitting) * lamb_nm_2pi, \
            'rms_photnoise': np.sqrt(sigma2_photnoise) * lamb_nm_2pi, \
            'rms_readnoise': np.sqrt(sigma2_readnoise) * lamb_nm_2pi, \
            'rms_anisop': np.sqrt(sigma2_anisop) * lamb_nm_2pi,'nphot':nphot,'Texp_ms':texp_ms}


    def compute_ee(self, psf, pixscale, fcore=1.5, spaxel=0):
        """
        :param strehl:
        :param seeing:
        :param filter:i
        :param teldiam:
        :param fcore:
        :return:
        """

        raper_ref = fcore * psf["FWHM_core"]

        # raper_core=raper_ref/sigma_core

        raper2_ref = raper_ref * raper_ref
        sigma2_halo = psf["sigma_halo"] * psf["sigma_halo"]
        sigma2_core = psf["sigma_core"] * psf["sigma_core"]
        print "========================="
        print raper2_ref.unit
        print sigma2_halo.unit
        print sigma2_core.unit
        print "========================="

        #psf = compute_psf(strehl, sigma_core, sigma_halo)

        A_core = psf['Amp_core']
        A_halo = psf['Amp_halo']
        fcore = psf['rat_core2halo']

        ee = (A_halo * sigma2_halo * (1 - np.exp(-raper2_ref.value / 2. / sigma2_halo.value)) + \
              A_core * sigma2_core * (1 - np.exp(-raper2_ref.value / 2. / sigma2_core.value))) / \
             (A_halo * sigma2_halo + A_core * sigma2_core)

        ee1 = 1. - 1. / (sigma2_core * fcore + sigma2_halo) * (sigma2_halo * np.exp(-raper2_ref.value / 2. / sigma2_halo.value) + \
                                                               fcore * sigma2_core * np.exp(-raper2_ref / 2. / sigma2_core))
        area_per_pixel = pixscale*pixscale
        ## check if spaxel is set, then assume a spaxel as 2x1 pixel, according to FRIDA specs.
        if (spaxel >0):
            area_per_pixel = area_per_pixel /2.

        area_apert=pi*raper_ref**2 ## in arcsec
        npix_apert = area_apert / area_per_pixel

        return {'EE': ee,'EE1': ee1, 'Radius':raper_ref, 'Npix':npix_apert, 'Area':area_apert,\
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


    def compute_psf(self,lamb_mic,strehl):
        """

        :param strehl: Strehl ratio
        :param sigma_core: Width of the core Gaussian (diffraction limit core)
        :param sigma_halo: Width of the halo Gaussian (seeing)
        :return:
        """
        fwhm_core = 0.0242 * lamb_mic * (10.4 / self.diam_telescope)
        ## seeing is assumed to be measured at 5000AA = 0.5 microns
        fwhm_halo = compute_seeing_lambda(lamb_mic,self.seeing_input,self.lambda_seeing_input)
        sigma_halo = fwhm_halo / 2.35
        sigma_core = fwhm_core / 2.35

        rat2_sigma = (sigma_halo / sigma_core) ** 2


        fcore = (strehl * rat2_sigma - (1 * (rat2_sigma.unit))) / (1 - strehl)


        amp_core_norm = 1. / 2 / pi / sigma_core ** 2 / (1 + rat2_sigma / fcore)

        amp_halo_norm = 1. / 2 / pi / sigma_core ** 2 / (fcore + rat2_sigma)

        return {"Amp_core": amp_core_norm, "Amp_halo": amp_halo_norm, "rat_core2halo": fcore,\
                "sigma_core":sigma_core,"sigma_halo":sigma_halo,"FWHM_core":fwhm_core,"FWHM_halo":fwhm_halo}

    def set_aperture_circular(self,pixscale,radius,spaxel=0):

        area_pixel = pixscale*pixscale
        if (spaxel >0):
            area_pixel = area_pixel /2.
        area_apert=pi*radius**2 ## in arcsec
        npix_apert = area_apert / area_pixel

        return {'radius':radius,'area_aperture':area_aperture,'Npix':npix_apert,'area_pixel':area_pixel}
