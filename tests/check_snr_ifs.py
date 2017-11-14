# load needed modules
from frida.set_config import *
from frida.calculator import SkyConditions,GuideStar,TargetInfo,Calculator_Spect
from frida.gtcao import GTC_AO
import astropy.units as u

from frida.compute_flux import *
import matplotlib.pyplot as plt


# sky conditions - provide by Get_SkyConditions(request)
sky_conditions={'seeing':0.9,'lambda_seeing':0.5 * u.micron,'airmass':1.2,'pwv':2.5}

# Telescope & instrument setup
data_path = settings.INCLUDES 

telescope =  'GTC'
telescope_params = Telescope(telescope)

instid = 'FRIDA'
scale = 'fine_scale'
instru_static = Instrument_static(scale=scale,instid=instid,path=data_path)


## call GTC_AO to compute Strehl ratio and Encircled Energy
gs_mag = 11.
gs_sep = 0.
guide_star = GuideStar(gs_mag,gs_sep)

aocor = GTC_AO(sky_conditions,guide_star,telescope_params.aperture)

grating_name='H Med'
cent_wave=1.65 * u.micron

dir_grating_list = settings.INCLUDES
dir_grating_info = settings.GRATINGS
dir_grating_skycalc = settings.SKYCALC_GRATINGS
grating_info = Grating(grating_name,path_list=dir_grating_list,path_gratings=dir_grating_info,path_skycalc=dir_grating_skycalc)

grating_info

print ("Dispersion ",grating_info.delt_wave)
print ("Sky rad:",grating_info.skyrad['value'][100:105])
#print ("Disp,lambda_center,bandwidth=",grating_info.dlambda,grating_info.lcenter,grating_width)
wave_array = grating_info.wave_array()

grating_effic = grating_info.interp_spect('GratingEffic',wave_array)
sky_rad = grating_info.interp_spect('SkyRad',wave_array)
print("sky_rad ",sky_rad[1000:1010])
atrans = grating_info.interp_spect('AtmosTrans',wave_array)

static_response = instru_static.compute_static_response(wave_array)
detector = instru_static.detector

throughput = static_response["collimator"] * static_response['camera'] * static_response['qe']
print ("Resp Coll ", static_response["collimator"][1000:1015])
print ("SPRF camera ", static_response["camera"][1000:1015])
print ("Detector QE ", static_response["qe"][1000:1015])
print ("Throughput ", throughput[1000:1015])

plt.plot(wave_array.to("angstrom"),static_response['qe'])
plt.plot(wave_array.to("angstrom"),static_response['collimator'])
plt.plot(wave_array.to("angstrom"),static_response['camera'])
plt.plot(wave_array.to("angstrom"),throughput)
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.title('Static response')
plt.savefig('tests/Static_response.png')
plt.close()

plt.plot(wave_array.to('angstrom'),grating_effic)
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
#plt.ylabel(r'$\rm{\ L_{\lambda} [{photons\ s^{-1} }]}$',size=20)
#plt.ylabel('r'+photunit.unit.to_string('latex'),size=20)
plt.savefig('tests/Grating_effic.png')
plt.close()


plt.plot(wave_array.to('angstrom'),atrans)
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
#plt.ylabel(r'$\rm{\ L_{\lambda} [{photons\ s^{-1} }]}$',size=20)
#plt.ylabel('r'+photunit.unit.to_string('latex'),size=20)
plt.savefig('tests/Atrans.png')
plt.close()

plt.plot(wave_array.to('angstrom'),sky_rad)
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
#plt.ylabel(r'$\rm{\ L_{\lambda} [{photons\ s^{-1} }]}$',size=20)
#plt.ylabel('r'+photunit.unit.to_string('latex'),size=20)
plt.savefig('tests/Sky_radiance.png')
plt.close()


## target info
print("Set target parameters")
mag_target = 18.
band = 'H'
mag_system='Vega'
sed=('black_body',3000.*u.K)
phot_zp = define_photzp()

print(phot_zp['J'])
print(phot_zp['J']['bwidth'])
# sed=('power_law',0.)
target_info = TargetInfo(mag_target,band,mag_system,sed)

flux_scale=target_info.flux_scale
print ("flux scale=",flux_scale)
photons_obj=target_info.photons_sed(wave_array)

photunit=u.photon/(u.s * u.cm**2 * u.angstrom)
print ("photons_obj=",photons_obj[1000:1010].to(photunit))
plt.plot(wave_array.to('angstrom'),photons_obj.to(photunit))
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
#plt.ylabel(r'$\rm{\ L_{\lambda} [{photons\ s^{-1} }]}$',size=20)
plt.ylabel('r'+photunit.to_string('latex'),size=20)
plt.savefig('tests/photons.png')
plt.close()

## # compute the effective photon rate from target and sky as measured by the detector
## (done within compute_photonrate_spect)
phi_obj_total = photons_obj * atrans * telescope_params.area * telescope_params.reflectivity \
                * grating_effic * throughput

electron_per_s_unit=u.electron/(u.s * u.angstrom)
exptime = 300 * u.s
print("Photons_detected_from_object in exptime ",exptime)
print("\t",phi_obj_total[1000:1010].to(electron_per_s_unit)*exptime)


## compute sky emission and atmospheric absorption  (done within compute_sky_photonrate_spect)
area_reference = np.pi*(50 * u.marcsec)**2
phi_sky_apert = sky_rad * telescope_params.area * telescope_params.reflectivity \
                * grating_effic * throughput * area_reference

print("Photons_detected_from_sky in exptime",exptime," reference area ",area_reference)
print("\t",phi_sky_apert[1000:1010].to(electron_per_s_unit)*exptime)


## compute S/N
dit = 10. * u.second
nexp = np.ceil((exptime / dit).value)
texp = nexp * dit
# print ("nexp ",nexp[0],nexp[-1])
aperture = {"EE":0.5,"Npix":27}
phi_obj_apert = phi_obj_total * aperture["EE"]
#noise = np.sqrt(texp * (phi_obj_apert + phi_sky_apert + detector['darkc'] * aperture['Npix']) + \
#                nexp * aperture['Npix'] * detector['ron'] ** 2)
print("detector_dark=",detector['darkc'])
print("detector_ron=",detector['ron'])
dwave = grating_info.delt_wave
shot_noise2 = (texp * dwave * (phi_obj_apert + phi_sky_apert).to(electron_per_s_unit)).value
dark_noise2 = (texp * detector['darkc'] * aperture['Npix']).value
read_noise2 = (nexp * aperture['Npix'] * detector['ron'] ** 2).value
print("shot noise=",np.sqrt(shot_noise2))
print("dark noise=",np.sqrt(dark_noise2))
print("read noise=",np.sqrt(read_noise2))
noise = np.sqrt(shot_noise2 + dark_noise2 + read_noise2)
#noise = np.sqrt(texp * (phi_obj_apert + phi_sky_apert ))
signal = texp * dwave * phi_obj_apert

print("Signal ",signal[1000:1010].to('electron'))
print("Noise ",noise[1000:1010])
snr = signal.to("electron").value / noise
print("SNR ",snr[1000:1010])

plt.plot(wave_array.to('angstrom'),snr)
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('S/N',size=20)
plt.savefig('tests/snratio.png')
plt.close()


###len(wave_array)
strehl= aocor.compute_strehl(wave_array[1000:1010],sky_conditions['airmass'])
print("Strehl ",strehl['StrehlR'])
psf = aocor.compute_psf(wave_array[1000:1010],strehl['StrehlR'])
## the aperture in the IFS case is box-like, one side corresponds to the "slit-width" and the other
## side is along the radial profile
fcore = 1.5 # radius of aperture as a factor of FWHM of the PSF core
pixscale = instru_static.pixscale # in arcseconds
ee_aperture=aocor.compute_ee(psf,pixscale,fcore=fcore,spaxel=1)
print("ee_aperture: ",ee_aperture)
