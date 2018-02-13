import numpy as np
import math

def compute_strehl(filter,seeing,mag_gs,dist_gs,airmass):
 print "AQUI"

 lamb_ref = 0.5
 r0_ref_cm = 20.*lamb_ref/seeing 
 diam_telescope = 10.4 # in meters
 theta_isoplan_ref = 2.0 # arcsec at 0.5 micron

 ## integration time is computed
 mag_gs_tmin = 13.5 
 c2_texp = 0.9
 if (mag_gs < mag_gs_tmin):
       texp_ms = 2.
 else:
       texp_ms = 2.+(mag_gs-mag_gs_tmin)**2*c2_texp        
 ## Parameters p0, p1 & p2 are extracted from table 1.3 RPT/OPTI/0252-R (GTCAO System Error Budgets) 3.A
  
 if (filter == 'H'):
      sr0_fit = 0.581
      ind_fit = 0.1637
      lamb_mic = 1.65
      noise_p0 = 0.2852
      noise_p1 = 8.8
      noise_p2 = 32. 
 elif (filter == 'J'):
      sr0_fit = 0.389
      ind_fit = 0.1644
      lamb_mic = 1.25
      rms_r0_20cm = 119
      noise_p0 = 0.4216
      noise_p1 = 16.5
      noise_p2 = 24. 
 elif (filter == 'K'):
      sr0_fit = 0.729
      ind_fit = 0.1468
      lamb_mic = 2.2
      noise_p0 = 0.1978
      noise_p1 = 4.8
      noise_p2 = 25. 

 r0_cm = r0_ref_cm * (lamb_mic/lamb_ref)**(6./5.)
 lamb_nm = lamb_mic*1.e3
 
 ## Scintillation
 h_cm = 5*1e5 
 lamb_cm = lamb_mic * 1.e-4
 sigma2_scint = 0.288*(np.sqrt(lamb_cm*h_cm)/r0_cm)**(5./3.)

 ## Static
 sigma_static_nm = 152.
 sigma2_static = (sigma_static_nm *2. * math.pi / lamb_nm)**2 
 

 ## bandwidth high-order
 f_G= 31.52
 f_s= 1./texp_ms*1.e3
 sigma2_ho = (f_G/f_s)**(5./3.)
 
 ## low-order tip-tilt or bandwidth low order (it uses a different expression)
 sigma_tt_bw_nrad = 11.3 ## due to bandwidth (specified in 1E-9 rad)
 sigma_tt_wndsh_nrad = 31.8 ## due to wind shake
 sigma_tt_nrad = sigma_tt_bw_nrad + sigma_tt_wndsh_nrad
 sigma_tt_rad = sigma_tt_nrad * 1.e-9 
 lamb_over_D = lamb_nm/diam_telescope
 strehl_tt= 1./(1+math.pi*math.pi/2.*(sigma_tt_nrad/lamb_over_D)**2)
 #f_T = 13.3
 #rms_tt_nm_bright = 110.   # reference value at 0.5 seeing y texp=2msec
 #f_s_ref = 1./2.e-3  ##  1/(texp_ref=2ms)
 #rms_tt_nm  = rms_tt_nm_bright * (f_s_ref/f_s)
 #sigma2_tt = (rms_tt_nm/lamb_nm*2*math.pi)**2  
 
 ## fitting error
 rms_r0_20cm = 119  # in nm 
 sigma2_fit= (rms_r0_20cm / lamb_nm * 2* math.pi) ** 2 * (20/r0_ref_cm)**(5./3.)  
   
 ## Photon noise error 
 effic_wfs = 0.27
 rad_subaper = 10.4/20./2.  ## diameter of primary divided into 20 subapertures 
 area_subaper = (rad_subaper)**2.*math.pi  ## area of subaperture in m^2
 nphot0 = 3.0216e7  ## photons/msec/m2/subaperture mag_R=0.
 nphot = nphot0 * 10**(-0.4*mag_gs) * texp_ms  * rad_subaper *effic_wfs
 sigma2_photnoise = noise_p1 / nphot 
 
 sigma2_readnoise = noise_p2 /nphot/nphot
 
 ## error due to separation of guide star relative to the pointing
 ##sigma_anisop = (dist_gs/theta_isoplan_ref)**(5./3)*(lamb_ref/lamb_mic)**2 
 sigma2_anisop = (dist_gs/theta_isoplan_ref)**(5./3.)*(lamb_ref/lamb_mic)**2

 ## adding all contributions computed before
 sigma2_total = sigma2_static + sigma2_scint + sigma2_ho + sigma2_fit + \
    sigma2_photnoise + sigma2_readnoise + sigma2_anisop
 
 ## includes the variation with airmass, which affects all sources of error 
 sigma2_total = sigma2_total * airmass
 
 strehl =np.exp(-sigma2_total)*strehl_tt
 lamb_nm_2pi = lamb_nm/2/math.pi
 
 ## as output produces Strehl ratio, RMS wavefront error components in nanometers 
 return {'StrehlR':strehl,'rms_scint':np.sqrt(sigma2_scint)*lamb_nm_2pi,'rms_stat':np.sqrt(sigma2_static)*lamb_nm_2pi,\
   'rms_bw_ho':np.sqrt(sigma2_ho)*lamb_nm_2pi,'strehl_tt':strehl_tt,\
   'rms_fitting':np.sqrt(sigma2_fit)*lamb_nm_2pi,\
   'rms_photnoise':np.sqrt(sigma2_photnoise)*lamb_nm_2pi,\
   'rms_readnoise':np.sqrt(sigma2_readnoise)*lamb_nm_2pi,\
   'rms_anisop':np.sqrt(sigma2_anisop)*lamb_nm_2pi}
 
def compute_ee(strehl,seeing,filter,teldiam=10.4,faper=1.5):

 if (filter == 'H'):
      lamb_mic = 1.65
 elif (filter == 'J'):
      lamb_mic = 1.25
 elif (filter == 'K'):
      lamb_mic = 2.2
      
 fwhm_core=0.0242*lamb_mic*(10.4/teldiam)
 ## seeing is assumed to be measured at 5000AA = 0.5 microns 
 fwhm_halo=seeing * (0.5/lamb_mic)**(1/5.)
 
 raper_ref=faper*fwhm_core

 sigma_halo=fwhm_halo/2.35
 sigma_core=fwhm_core/2.35
 #raper_core=raper_ref/sigma_core

 raper2_ref=raper_ref*raper_ref
 sigma2_halo=sigma_halo*sigma_halo
 sigma2_core=sigma_core*sigma_core
 
 psf=compute_psf(strehl,sigma_core,sigma_halo)

 A_core=psf['Amp_core']
 A_halo=psf['Amp_halo']
 fcore=psf['rat_core2halo']
 

 ee =(A_halo*sigma2_halo*(1-np.exp(-raper2_ref/2./sigma2_halo))+\
    A_core*sigma2_core*(1-np.exp(-raper2_ref/2./sigma2_core)))/\
    (A_halo*sigma2_halo+A_core*sigma2_core)  
 
 ee1 = 1.- 1./(sigma2_core*fcore+sigma2_halo)*(sigma2_halo*np.exp(-raper2_ref/2./sigma2_halo)+\
         fcore*sigma2_core*np.exp(-raper2_ref/2./sigma2_core))       
      
 return {'EE':ee,'FWHM_halo':fwhm_halo,'FWHM_core':fwhm_core,'Amp_halo':A_halo,'Amp_core':A_core,'EE1':ee1}  

def compute_psf(strehl,sigma_core,sigma_halo):
  
  rat2_sigma = (sigma_halo/sigma_core)**2 
  
  fcore = (strehl*rat2_sigma-1)/(1-strehl)
  
  amp_core_norm = 1./2/math.pi/sigma_core**2/(1+rat2_sigma/fcore)
  
  amp_halo_norm = 1./2/math.pi/sigma_core**2/(fcore+rat2_sigma)
  
  return {"Amp_core":amp_core_norm,"Amp_halo":amp_halo_norm,"rat_core2halo":fcore} 
