import numpy as np
import matplotlib.pyplot as plt
from frida.set_config import *
from frida.gtcao import *
from frida.aux_functions import *
from frida.calculator import GuideStar
import scipy.constants as pconst
import astropy.units as u
from frida.calculator import SkyConditions, GuideStar, TargetInfo, Calculator_Image, Calculator_Spect


# telescope 
telescope = "GTC"
telescope_params = Telescope(telescope)

# Atmospheric conditions - fried parameter
#r0=10   
fwhm_seeing=0.9
lambda_seeing=0.55
airmass=1.2

## guide star characteristics
mag_gs=13.
sep_gs=0. * u.arcsec

guide_star = GuideStar(mag_gs,sep_gs)

frida_static = Instrument_static()

print('instrument_detector:',frida_static.detector)
#pixscale = 0.01 * u.arcsecond # in arcseconds
pixscale = frida_static.pixscale

print('pixscale',pixscale)
print('PSF_model ',settings.PSF_MODEL)

sky_conditions = {'seeing': 0.9 * u.arcsec, 'wave_seeing': 0.5 * u.micron, \
        'airmass': 1.2, 'pwv': 2.5 * u.mm}

aocor = GTC_AO(sky_conditions,guide_star,telescope_params.aperture)

wave_input = np.array([5000.,11000,21000.]) * u.Angstrom

r0_wave = compute_r0(wave_input,sky_conditions['seeing'],\
                     wave_ref=sky_conditions['wave_seeing'])

print('r0_wave=',r0_wave.to('cm'))

#wave_in = [11500.,16500.,22000.] * u.angstrom
wave_in = 22000. * u.angstrom

strehl=aocor.compute_strehl(wave_in,airmass)


mag_target = 18.
mag_system = 'Vega'
band ='K'
temperature = 3000.
sed = ('black_body',temperature)
extended=False
label_source_morpho = 'Point source'
waveshift = ('radial_velocity', 1. * u.km / u.s )
label_energy_type = 'Black Body, T=%s K' % temperature
target_info = TargetInfo(mag_target,band,mag_system,sed,extended=extended)
#target_info.error = error
#target_info.messages = debug_values
#target_info.error_messages = error_values
target_info.energy_type = label_energy_type
target_info.source_type= label_source_morpho

## 
calc_filter = 'Ks'
calc_scale = 'fine_scale'
a=Calculator_Image(target_info,calc_filter,scale=calc_scale)
unit_e_sec = u.electron/u.second

print("a.phi_obj_total: ",a.phi_obj_total.to(unit_e_sec))
print("a.Sky emission mag: ",a.skymag) 
print("a.Sky emission photons: ",a.skyemission_photons) 
print("a.Sky emission photons/arcsec^2/s: ",a.phi_sky_sqarc) 
print("a_pixscale: ",a.pixscale) 
print("a.Sky emission photons/pixel/s: ",a.phi_sky_sqarc * a.pixscale**2)

ron = frida_static.detector['ron']
darkc = frida_static.detector['darkc']
dit = 20. * u.second


wave_eff = wave_input[2]
psf_2gauss = aocor.compute_psf('2-Gaussians',wave_eff,strehl['StrehlR'])
print("PSF-2Gauss: ",psf_2gauss)
psf_airy = aocor.compute_psf('Airy+Gaussian',wave_eff,strehl['StrehlR'])
print("PSF-Airy: ",psf_airy)

Nx=200
Ny=200
psf2d_2gauss = buildim_psf2d_2gauss(psf_2gauss,pixscale,Nx=Nx,Ny=Ny) 
psf2d_airy = buildim_psf2d_AiryGauss(psf_airy,pixscale,Nx=Nx,Ny=Ny) 

im2show = (psf2d_airy* (a.phi_obj_total*dit).to(u.electron)).value
print("im2show [100,100]",im2show[100,100])
vmin = np.median(im2show)*0.5
vmax = np.median(im2show)*100.5
plt.imshow(im2show,vmax=vmax,vmin=vmin)
plt.colorbar()
plt.show()

atmphere = Atmosphere()
atmos_skyrad = interpol2newwave(atmphere.skyrad_photons,\
                atmphere.skyrad_wave,wave_eff)
area_pixel = a.pixscale*a.pixscale
sky_pixel = a.phi_sky_sqarc * area_pixel 
print('sky_pixel:',sky_pixel) 

psf_object =psf2d_airy * a.phi_obj_total.to(unit_e_sec) * dit  
print("psf_object [100,100]",psf_object[100,100])
psf_object_sky =psf_object + sky_pixel*dit
print("psf_object_sky [100,100]",psf_object_sky[100,100])


im2show = psf_object_sky.value    
print("im2show [100,100]",im2show[100,100])
vmin = np.median(im2show)*0.5
vmax = np.median(im2show)*1.5
print("median,vmin,vmax:",np.median(im2show),np.max(im2show),vmin,vmax)
plt.imshow(im2show,vmax=vmax,vmin=vmin)
plt.colorbar()
plt.show()
    
# now compute the noise
im_obj = psf2d_airy * a.phi_obj_total.to(unit_e_sec) 
noise_2d= build_im_signal_with_noise(im_obj,sky_pixel,\
                ron,darkc,dit,Nexp=1,SubtractBck=True)

Nexp =1
noise_psf_object_1d = np.random.poisson(lam=(psf_object_sky.value).reshape(Nx*Ny),size=(1,Nx*Ny))
noise_darkc_2d = (np.random.poisson(lam=(dit*darkc).value,size=(Nx*Ny))).reshape(Nx,Ny)
noise_ron_2d = (np.random.normal(0.,np.sqrt(2*Nexp)*ron.value,size=(Nx*Ny))).reshape(Nx,Ny)
print("RON",dit)
print("RON*Nexp",np.sqrt(2*Nexp)*ron.value)

print("RON min,max,median",np.min(noise_ron_2d),np.max(noise_ron_2d),np.median(noise_ron_2d))
print("Darkc min,max,median",np.min(noise_darkc_2d),np.max(noise_darkc_2d),np.median(noise_darkc_2d))

noise_psf_object_2d = noise_psf_object_1d.reshape(Nx,Ny) 
print("shape noise_psf_object_2d=",noise_psf_object_2d.shape)
print("median noise_psf_object_2d=",np.median(noise_psf_object_2d))

#point_with_noise_2d=build_im_signal_with_noise(psf_object,sky_pixel,\
#                ron,darkc,dit,Nexp=1,SubtractBck=True)
im2show = noise_psf_object_2d +noise_darkc_2d + noise_ron_2d - (sky_pixel*dit).value
print("im2show [100,100]",im2show[100,100])
vmin = np.median(im2show)*0.8
vmax = np.min([np.max(im2show)*0.5,np.median(im2show)*15.])
print("median,vmin,vmax:",np.median(im2show),np.max(im2show),vmin,vmax)
plt.imshow(im2show,vmax=vmax,vmin=vmin)
plt.colorbar()
plt.show()

'''
fcore = 1.2
aperture=aocor.compute_ee(psf_airy,pixscale,fcore=fcore)
print("fcore ",fcore,"   Aperture:",aperture)

fcore = 10.5
aperture=aocor.compute_ee(psf_airy,pixscale,fcore=fcore)
print("fcore ",fcore,"   Aperture:",aperture)

## instrument characteristics
filter='J'
#param_filter = get_filter(filter)
#lambda_obs = param_filter['lamb_eff']

## compute strehl and ee
#strehl = compute_strehl(filter,fwhm_seeing,mag_gs,sep_gs,airmass)
#ee_point = compute_ee(strehl['StrehlR'],fwhm_seeing,filter)
#radaper = 1.5*ee_point['FWHM_core']
#area_ref = pconst.pi*radaper*radaper  # arcsec^2
#npix =  area_ref / pixscale / pixscale
'''

 
