from math import exp, pow, sqrt
import numpy as np
import random, os, glob

import astropy.units as u
import astropy.constants as const
from astropy.modeling.blackbody import blackbody_lambda
#from astropy.analytic_functions import blackbody_lambda
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

#import frida.gtcao
#from frida.gtcao import *
#import frida.compute_flux
from frida.compute_flux import *
#import frida.set_config
from frida.set_config import *
from frida.aux_functions import *

phot_zp = define_photzp()


class GuideStar:
	""" Class containing information about Guide Star characteristics (magnitude and distance to opt. axis)
    :param gs_mag: Magnitude of guide Star 
    :type gs_mag: float. 
    :param gs_dist: Distance of guide Star to optical axis. 
    :type gs_dist: float 
    :returns: class self. -- .Magnitude, .Separation 
    :raises: 
	"""
	def __init__ (self, gs_mag=8, gs_dis=0.*u.arcsec):
		self.Magnitude = gs_mag
		self.Separation = gs_dis


class TargetInfo:
	""" Class containing information about Science Target
    :param mag_target: Magnitude of guide Star 
    :type mag_target: float. 
    :param band: Band in which the brightness is specified. 
    :type band: string 
    :param mag_system: Magnitude system (Vega or AB). 
    :type mag_system: string 
    :param sed: Type of Spectral Energy Distribution 
    :type sed: string 
    :returns: class self. -- .Magnitude, .Separation 
    :raises: 
	"""
	def __init__ (self,mag_target,band,mag_system,sed,waveshift=('radial_velocity', 0.), wave_ini=8500*u.angstrom,wave_end=25000.*u.angstrom,dwave=2.*u.angstrom):
		self.Magnitude = mag_target
		self.Band = band
		self.MagSystem = mag_system
		self.SED = sed
		self.lambda_band= phot_zp[band]['efflamb']   # calculamos el flujo de referencia

		## waveshift es una tupla el primer elemento es una cadena indicando \
		## la seleccion "radial_velocity" o "redshift", el segundo elemento es 
		## el valor de la velocidad, que debe ser en km/s o el redshift. 
		## Por ejemplo waveshift =('radial_velocity',350.) o ('redshift',0.084)
		redshift = float(waveshift[1])
		if (waveshift[0] == "radial_velocity"):
			redshift = redshift*u.Unit('km/s')/const.c.to('km/s')
		lambda_band_rest = self.lambda_band/(1.+redshift)
			
        
		##REFACTOR --- define functions 
		if (sed[0] == "black_body"):
			temp_bb = float(sed[1])*u.K
			# calculamos el flujo en la banda de referencia y el factor de escala
			#input_normflux=bbody(lambda_band_rest,temp_bb)
			omega=1.e-16*u.sr
			input_normflux=blackbody_lambda(lambda_band_rest,temp_bb)*omega
			flux_scale=scale_to_vega(band,mag_target,input_normflux)
			wave_array_rest = np.arange(wave_ini.to('angstrom').value,\
			   wave_end.to('angstrom').value,dwave.to('angstrom').value)\
               * u.angstrom
			print(' normflux=',input_normflux,' value=',input_normflux.value)  
			print(' flux_scale=',flux_scale,' value=',flux_scale.value,float(flux_scale))   
			#sed_flambda = flux_scale * bbody(wave_array_rest,temp_bb)
			sed_wave = wave_array_rest*(1+redshift)
			sed_flambda = flux_scale * blackbody_lambda(wave_array_rest,temp_bb)*omega
		elif (sed[0] == "power_law"):
			index_pwl = sed[1]
			input_normflux=powl(self.lambda_band,index_pwl)
			flux_scale=scale_to_vega(band,mag_target,input_normflux)
			print(' normflux=',input_normflux,' value=',input_normflux.value)  
			print(' flux_scale=',flux_scale,' value=',flux_scale.value,float(flux_scale))   
			wave_array_rest = np.arange(wave_ini.to('angstrom').value,\
			   wave_end.to('angstrom').value,dwave.to('angstrom').value)\
               * u.angstrom
			sed_wave = wave_array_rest*(1+redshift)
			sed_flambda = flux_scale * powl(sed_wave,index_pwl)
		elif (sed[0] == "pickles"):
		    template = sed[1]
		    #sed_pickles=read_sed_pickles_files(template,\
		    #      path_sed=os.path.join(settings.SED_LIBRARY,'pickles'))
		    sed_pickles=read_sed_pickles_files(template,\
		          path_sed=os.path.join(settings.SED_LIBRARY,'pickles'))
		    print("sed_pickles_wave=",sed_pickles['wave'])
		    print("sed_pickles_flambda=",sed_pickles['flambda'])
		    print("lambda_band_rest=",lambda_band_rest)
		    lambda_band_rest_sameunit = lambda_band_rest.to(sed_pickles['wave'].unit).value
		    input_normflux=np.interp(lambda_band_rest_sameunit,sed_pickles['wave'].value,\
		          sed_pickles['flambda'].value) * sed_pickles['flambda'].unit
		    print("input_normflux=",input_normflux)
		    flux_scale = scale_to_vega(band,mag_target,input_normflux) 
		    print("flux_scale=",flux_scale)
		    sed_flambda = flux_scale * sed_pickles['flambda']
		    sed_wave = sed_pickles['wave']*(1+redshift)
		elif (sed[0] == "nonstellar"):
		    template = sed[1]
		    sed_nonstellar=read_sed_nonstellar_files(template,\
		          path_sed=os.path.join(settings.SED_LIBRARY,'nonstellar'))
		    print("sed_nonstellar_wave=",sed_nonstellar['wave'])
		    print("sed_nonstellar_flambda=",sed_nonstellar['flambda'])
		    print("lambda_band_rest=",lambda_band_rest)
		    lambda_band_rest_sameunit = lambda_band_rest.to(sed_nonstellar['wave'].unit).value
		    input_normflux=np.interp(lambda_band_rest_sameunit,sed_nonstellar['wave'].value,\
		          sed_nonstellar['flambda'].value) * sed_nonstellar['flambda'].unit
		    flux_scale = scale_to_vega(band,mag_target,input_normflux) 
		    sed_flambda = flux_scale * sed_nonstellar['flambda']
		    sed_wave = sed_nonstellar['wave']*(1+redshift)
		## END Refactor

		print(" min - max sed_wave",min(sed_wave),max(sed_wave))
		mask = (sed_wave >= wave_ini) & (sed_wave <= wave_end) 
		self.flux_scale = flux_scale
		self.sed_wave = sed_wave[mask]
		self.sed_flambda = sed_flambda[mask]

	def debug(self):
		result = []
		result.append("=========== DEBUG TARGET INFO ==============")
		result.append("Magnitude: %.2lf" % self.Magnitude)
		result.append("Band: %s" % self.Band)
		result.append("Mag System: %s " % str(self.MagSystem))
		result.append("SED: %s" % str(self.SED))
		result.append("Lambda Band: %s" % str(self.lambda_band))   # calculamos el flujo de referencia)
		result.append("Flux Scale: %s" % str(self.flux_scale))
		result.append("Sed Wave: %s" % str(self.sed_wave))
		result.append("Sed Flambda: %s" % str(self.sed_flambda))
		result.append("===========================================")
		print ('\n'.join(result))
		return result


	def flambda_wave(self,wave):
		"""
		Compute f_lambda for the selected SED and a given wave array
		:param wave:
		:return:
		"""

		#if (self.SED[0] == "black_body"):
		#	flambda = self.flux_scale * bbody(wave,self.SED[1])
		#elif (self.SED[0] == "power_law"):
		#	flambda = self.flux_scale * powl(wave,self.SED[1])
        ## interpolate 
		flambda_unit = self.sed_flambda.unit
		if hasattr(wave,'unit'):
			wave_sameunit = wave.to(self.sed_wave.unit).value
		else:
			wave = wave * u.micron
			wave_sameunit = wave.to(self.sed_wave.unit).value
		flambda = np.interp(wave_sameunit,self.sed_wave.value,self.sed_flambda.value)
		return flambda * flambda_unit

	def photons_wave(self,wave):
		"""
		Compute f_lambda for the selected SED and a given wave array
		:param wave:
		:return:
		"""
		phot_wave = self.flambda_wave(wave)/energy_photon(wave) #const.h/const.c*wave
		return phot_wave


#zeropoint=  {
#					"I":[0.79, 7.91, 7.77],
#					"J":[1.25, 8.51, 8.14],
#					"H":[1.65, 8.94, 8.39],
#					"K":[2.16, 9.40, 8.64]
#				}

class SkyConditions:

	Seeing = None
	AirMass = None

	def __init__ (self, seeing = 0.75, airmass = 1.0):
		self.Seeing = seeing
		self.AirMass = airmass


### added by JAP

class Calculator_Image:

	energy_photon_1mic = 1.985e-19  ### Watt=[kg m^2 s^-2]

	obs_filter = None
	atmosphere = None

	def __init__(self,target_info,filter_name,atmospheric_cond=None,scale=None,\
              telescope_name='GTC',instrument_name='FRIDA'):
		"""
			Celestial source parameters in a dictionary
			target_info = {'Magnitude':19,'MagSystem':,'Band':'J','sed':sed} ;
			 sed is a tuple like ('PowerLaw',index_powerlaw), ('Blackbdoy',temperature_bb)
		"""
		self.debug_values = {}
		"""
		sed = target_info.SED
		temp_bb = sed[1]
		input_band = target_info.Band   # definimos como banda de referencia
		input_lambda= phot_zp[input_band]["bwidth"]   # calculamos el flujo de referencia
		mag_ref = float(target_info.Magnitude)
		input_normflux=bbody(input_lambda,temp_bb)  # calculamos el flujo en la banda de referencia
		flux_scale =scale_to_vega(input_band,mag_ref,input_normflux)
		flux_scale=target_info.flux_scale
		print ("flux scale",flux_scale)
		"""


		# telescope
		telescope = Telescope(telescope_name)
		self.area_tel = telescope.area
		self.refl_tel = telescope.reflectivity
		print ('Telescope is ',telescope.name)
		print ("  reflectivity is",self.refl_tel)


		#### filter
         ### wave,transmission are arrays with the bandpass
		self.obs_filter_name=filter_name
		obs_filter = Filter(filter_name,path_list=settings.INCLUDES,path_filters=settings.FILTERS)
		lambda_center = obs_filter.lambda_center()
		print("lambda_effective ",lambda_center)
		filter_width = obs_filter.width()
		print ("lambda_eff,bandwidth=",lambda_center,filter_width["width"],filter_width["cut_on"],filter_width["cut_off"])
		self.filter_transavg = obs_filter.avg_transmission()
		print ("filter_transmission_avg",obs_filter.avg_transmission())

		# Atmosphere
		atmphere = Atmosphere()

		## do in nm
		same_unit = obs_filter.wave_limits[0].unit
		dwave_img  = min([atmphere.atmtrans_delta.to(same_unit),atmphere.skyrad_delta.to(same_unit)])

		wave_img = np.arange(obs_filter.wave_limits[0].value,obs_filter.wave_limits[1].value,\
                       dwave_img.value)
		wave_img *= same_unit

		atmos_trans = interpol2newwave(atmphere.atmtrans_perone,\
                atmphere.atmtrans_wave,wave_img)

		atmos_skyrad = interpol2newwave(atmphere.skyrad_photons,\
                atmphere.skyrad_wave,wave_img)

		filter_trans = interpol2newwave(obs_filter.transmission,\
                obs_filter.wave,wave_img)

		self.atmostrans = atmos_trans
		self.skyemission_photons = atmos_skyrad
		self.skymag = atmphere.compute_skymag(filter_name)

		self.img_wave = wave_img
		self.img_filter_trans = filter_trans


		## instrument characteristics
		instrument=Instrument_static(scale,instid=instrument_name)
		#self.instrument = instrument
		static_response = instrument.compute_static_response(wave_img)
		self.static_response = static_response
		self.throughput = static_response["collimator"] * static_response['camera'] \
            * static_response['qe'] * telescope.reflectivity
		self.detector = instrument.detector
		self.pixscale = instrument.pixscale
        


		## compute the target f-lambda within the wavelength range determined by the filter transmission
		obj_flambda_spec = target_info.flambda_wave(wave_img)
		obj_photons_spec = target_info.photons_wave(wave_img)
		ii = int(len(wave_img)/2)
		print("QE ",static_response['qe'][ii-5:ii+5])
		print("Collimator ",static_response['collimator'][ii-5:ii+5])
		print("Camera ",static_response['camera'][ii-5:ii+5])
		print (" wave,obj_flambda,photons["+str(ii)+"]=",wave_img[ii],obj_flambda_spec[ii], obj_photons_spec[ii])

		'''

		## compute the flux from the target and the number of photons per unit time
		#self.obj_flambda = flux_scale * bbody(filter_info['lamb_eff'],temp_bb) # Units are
		'''
		self.obj_flambda_lambcenter = target_info.flambda_wave(lambda_center)
		print ("obj_flambda_lambeff:",self.obj_flambda_lambcenter)

		avg_response = self.get_average_response()
		print("avg_response[effic_static] ",avg_response["effic_static"])        
		print("avg_response ",avg_response)        
		phi_obj_total2 = self.compute_photonrate_simple(self.obj_flambda_lambcenter,\
                lambda_center,filter_width["width"],avg_response["atm_trans"]*\
                avg_response["effic_static"])

		self.flambda_wave = target_info.flambda_wave(wave_img)
		self.phi_obj_total = self.compute_target_photonrate_filter(target_info.flambda_wave(wave_img))

		units_phi = u.electron / u.second
		print("phi_obj_total [simple int-filter]",phi_obj_total2.to(units_phi),\
            self.phi_obj_total.to(units_phi))


		## compute the photon rate from the sky per arcsec^2
		##sky_flambda = 10.**(phot_zp[self.obs_filter_name][2]-0.4*self.atmosphere['mag_sqarc'])
		sky_flambda = self.get_skyflux_band(lambda_center)
		phi_sky_sqarc_simple = \
              self.compute_photonrate_simple(avg_response['sky_rad'],\
              lambda_center,filter_width["width"],avg_response["effic_static"],\
              is_flux=False)
		self.phi_sky_sqarc = self.compute_sky_photonrate_filter(atmos_skyrad)

		print("phi_sky/pix/sec (simple  int-filter)",phi_sky_sqarc_simple,\
			  self.phi_sky_sqarc)
        
		print("pixscale ",instrument.pixscale)
         
        	
	def get_average_response(self):
		'''
           compute quantities averaged over filter transmission 
           self.img_filter_trans
           self.img_wave
           self.static_response["collimator"]
		''' 
		coll_transavg = average_bandpass(self.static_response["collimator"],\
                                   self.img_filter_trans,self.img_wave)
		cam_transavg = average_bandpass(self.static_response["camera"],\
                                   self.img_filter_trans,self.img_wave)
		qe_avg = average_bandpass(self.static_response["qe"],\
                                   self.img_filter_trans,self.img_wave)
		effic_avg = self.refl_tel * coll_transavg * \
				  self.filter_transavg * cam_transavg * qe_avg
		atm_transavg = average_bandpass(self.atmostrans,\
                                   self.img_filter_trans,self.img_wave)
		sky_radavg = average_bandpass(self.skyemission_photons,\
                                   self.img_filter_trans,self.img_wave)        
		return {"collim":coll_transavg,"camera":cam_transavg,"qe":qe_avg,\
           "effic_static":effic_avg,"atm_trans":atm_transavg,"sky_rad":sky_radavg}

	def get_static_response(self):
		return (self.static_response, self.throughput)

	def compute_photonrate_simple(self,flux_filter,lambda_filter,width_filter,\
                    throughput,is_flux=True):
		"""
		compute the number of photons
		INPUTS:
		"""
		"""
		   Units: [flux_band] is expected in J s^-1 m^-2 micron^-1
				  [h] = J s
				  [lambda_band,width_band] = micron
				  [telarea] = m^2
				  [throughput] = No units
					Results as ===> s^-1 m^-1 micron -> x1.E-6 -> s^-1
		"""
		phot_tot = flux_filter*width_filter*self.area_tel*throughput  # units are photons/s
		if (is_flux): phot_tot /= energy_photon(lambda_filter)
		return phot_tot

	def compute_target_photonrate_filter(self,flambda,Atm_absorption=True):
		"""
		compute the number of photons integrated upon filter passband
		INPUTS:
		"""
		"""
		   Units: [flux_band] is expected in J s^-1 m^-2 micron^-1
				  [h] = J s
				  [lambda_band,width_band] = micron
				  [telarea] = m^2
				  [throughput] = No units
					Results as ===> ph^-1 s^-1 m^-1 micron -> x1.E-6 -> ph^-1 s^-1
		"""

		#a = self.integrand_obj(flambda) # include convolution with atmospheric and filter transmission
		#self.debug_values['a'] = a[0]
		#self.debug_values['flambda'] = flambda[0]
		photons_wave = flambda * self.img_filter_trans * \
            self.throughput / energy_photon(self.img_wave)
		if (Atm_absorption): photons_wave *= self.atmostrans 
         
		photons_through_filter = self.sum_photons_filter(photons_wave)
		#scipy.integrate.simps(integrand_obj.value,self.img_wave.value)
		#photons_through_filter = scipy.integrate.simps(integrand_obj.value,self.img_wave.value)
		#units_photons_trough_filter = integrand_obj.unit * self.img_wave.unit
		#return (photons_through_filter*self.area_tel*self.refl_tel).value * u.ph * u.s**-1
		return photons_through_filter*self.area_tel

	def compute_sky_photonrate_filter(self,photon_flux):
		"""
		compute the number of photons integrated upon filter passband
		INPUTS:
		"""
		"""
		   Units: [flux_band] is expected in J s^-1 m^-2 micron^-1
				  [h] = J s
				  [lambda_band,width_band] = micron
				  [telarea] = m^2
				  [throughput] = No units
					Results as ===> ph^-1 s^-1 m^-1 micron -> x1.E-6 -> ph^-1 s^-1
		"""

		#a = self.integrand_obj(flambda) # include convolution with atmospheric and filter transmission
		#self.debug_values['a'] = a[0]
		#self.debug_values['flambda'] = flambda[0]
		photons_wave = photon_flux * self.img_filter_trans * self.throughput 
         
		photons_through_filter = self.sum_photons_filter(photons_wave)
		#scipy.integrate.simps(integrand_obj.value,self.img_wave.value)
		#photons_through_filter = scipy.integrate.simps(integrand_obj.value,self.img_wave.value)
		#units_photons_trough_filter = integrand_obj.unit * self.img_wave.unit
		#return (photons_through_filter*self.area_tel*self.refl_tel).value * u.ph * u.s**-1
		return photons_through_filter*self.area_tel

	def sum_photons_filter(self,photons):
		photons_filter = scipy.integrate.simps(photons.value,self.img_wave.value)
		units_photons_trough_filter = photons.unit * self.img_wave.unit
		return photons_filter * units_photons_trough_filter


	def signal_noise_texp_img(self,Nexp_min,Nexp_max,Nexp_step,dit,aperture):
		"""
		Compute the S/N ratio as function of exposure time
		:param Nexp_min:
		:param Nexp_max:
		:param Nexp_step:
		:param dit:
		:param aperture:
		:return:
		"""
		nexp = np.arange(Nexp_min,Nexp_max,Nexp_step)
		texp = nexp * dit
		phi_obj_apert = self.phi_obj_total * aperture['EE']
		phi_sky_apert = self.phi_sky_sqarc * aperture['Area']

		unit_signal = u.electron

		noise2_obj = dit*phi_obj_apert
		noise2_sky = dit*phi_sky_apert
		noise2_read = aperture['Npix']*self.detector['ron']**2
		noise2_dark = dit*self.detector['darkc']
		print("===========================")
		print ("noise2_obj ",noise2_obj.to(unit_signal))
		print ("noise2_sky ",noise2_sky.to(unit_signal))
		print ("noise2_read ",noise2_read,noise2_read.value,noise2_read.unit)
		print ("noise2_dark ",noise2_dark,noise2_dark.value,noise2_dark.unit)
		print ("===========================")
		noise = np.sqrt(nexp*(noise2_obj.to(unit_signal).value+\
                        noise2_sky.to(unit_signal).value+noise2_read.value+\
                        noise2_dark.value))
        
		##noise = np.sqrt(texp*(phi_obj_apert+phi_sky_apert+self.detector['darkc']*aperture['Npix'])+ \ nexp*aperture['Npix']*self.detector['ron']**2)
		self.debug_values['phi_sky_sqarc'] =self.phi_sky_sqarc
		self.debug_values['phi_obj_total'] = self.phi_obj_total
		self.debug_values['phi_obj_aperture'] = phi_obj_apert.to(u.electron/u.second)
		self.debug_values['phi_sky_aperture'] = phi_sky_apert
		signal = texp * phi_obj_apert
		snr = signal.to(unit_signal).value / noise
		print("@signal_noise_texp_img - snr ",snr) 
		self.debug_values['signal'] = signal.to(unit_signal)
		#self.debug_values['noise'] = (noise * signal.unit).to(unit_signal)
		self.debug_values['noise'] = noise * unit_signal
		self.debug_values['snr'] = snr
		print("@signal_noise_texp_img - noise ",self.debug_values['noise'][10:15]) 
    
		return (texp,snr)

	def texp_signal_noise_img(self,required_sn,dit,aperture):
		"""
		Compute the exposure time needed to reach a S/N ratio
		:param Nexp_min:
		:param Nexp_max:
		:param Nexp_step:
		:param dit:
		:param aperture:
		:return:
		"""
		phi_obj_apert = self.phi_obj_total * aperture['EE']
		phi_sky_apert = self.phi_sky_sqarc * aperture['Area']

		noise2_obj = dit*phi_obj_apert
		noise2_sky = dit*phi_sky_apert
		noise2_read = aperture['Npix']*self.detector['ron']**2
		noise2_dark = dit*self.detector['darkc']

		unit_signal = u.electron

		print ("===========================")
		print (phi_obj_apert.unit)
		print (phi_sky_apert.unit)
		print (self.detector['darkc'],"***",aperture['Npix'])
		print (self.detector['ron']**2,"****",aperture['Npix'])
		print ("===========================")
		noise2 = noise2_obj.to(unit_signal).value+\
                        noise2_sky.to(unit_signal).value+noise2_read.value+\
                        noise2_dark.value
		signal2 = (phi_obj_apert * dit)**2
		print("signal2 ",signal2)

		nexp = round(required_sn**2 * noise2 / signal2.to(unit_signal*unit_signal).value)

		return nexp


	def get_skyflux_band(self,wave):
		from frida.aux_functions import interpolate
		wave_zp = []
		photzp_values = []
		matrix = []

#		for each_band in phot_zp:
#			wave_zp.append(phot_zp[each_band]['bwidth'])
#			photzp_values.append(phot_zp[each_band]['mAB_to_Vega'])
#
#		if hasattr(wave_zp[0],'unit'):
#			wave_zp = [item.value for item in wave_zp]
#		if hasattr(photzp_values[0],'unit'):
#			photzp_values = [photzp.value for photzp in photzp_values]
#
#		phot_zp_interp = interpolate(np.array(wave_zp),np.array(photzp_values),wave)
#
#		skyflux = 10.**(phot_zp_interp-0.4*self.skymag)

		for each_band in phot_zp:
			if hasattr(phot_zp[each_band]['bwidth'], 'unit'):
				current_wave = phot_zp[each_band]['bwidth'].value
			else:
				current_wave = phot_zp[each_band]['bwidth']

			if hasattr(phot_zp[each_band]['mAB_to_Vega'], 'unit'):
				current_mag = phot_zp[each_band]['mAB_to_Vega'].value
			else:
				current_mag = phot_zp[each_band]['mAB_to_Vega']

			print("## ", each_band, current_wave)

			matrix.append([current_wave, current_mag])
		matrix = np.sort(np.array(matrix), axis=0)
		phot_zp_interp = interpolate(matrix[:,0], matrix[:,1], wave)
		skyflux = 10.**(phot_zp_interp-0.4*self.skymag)


		return skyflux



class Calculator_Spect:

	energy_photon_1mic = 1.985e-19  ### Watt=[kg m^2 s^-2]

	obs_filter = None
	atmosphere = None

	def __init__(self,target_info,sky_conditions,spectrograph_setup,scale,Atmospheric_Cond=None,telescope=None,instrument=None, aocor=None):
		"""
			Celestial source parameters in a dictionary
			target_info = {'Magnitude':19,'MagSystem':,'Band':'J','sed':sed} ;
			sed is a tuple like ('PowerLaw',index_powerlaw), ('Blackbdoy',temperature_bb)
		"""
		# Setup the wavelength array from spectrograph_setup
		#### grating
		self.grating_info = Grating(spectrograph_setup['Grating'],cent_wave=spectrograph_setup['Central_wave'])
		self.telescope_params = telescope
		self.instrument = instrument
		self.aocor = aocor

	def get_static_response(self):
		wave_array = self.grating_info.wave_array()
		static_response = self.instrument.compute_static_response(wave_array)
		self.static_response = static_response
		self.throughput = static_response["collimator"] * static_response['camera'] * static_response['qe']
		return (self.static_response, self.throughput)

	def get_atrans(self):
		wave_array = self.grating_info.wave_array()
		atrans = self.grating_info.interp_spect('AtmosTrans',wave_array)
		return (atrans)

	def get_sky_rad(self):
		wave_array = self.grating_info.wave_array()
		sky_rad = self.grating_info.interp_spect('SkyRad',wave_array)
		return (sky_rad)

	def get_grating_effic(self):
		wave_array = self.grating_info.wave_array()
		grating_effic = self.grating_info.interp_spect('GratingEffic',wave_array)
		return (grating_effic)
		



	def compute_sky_photonrate_spect(self):
		"""
		compute the number of photons integrated upon filter passband
		INPUTS:
		   Units: [flux_band] is expected in J s^-1 m^-2 micron^-1
				  [h] = J s
				  [lambda_band,width_band] = micron
				  [telarea] = m^2
				  [throughput] = No units
					Results as ===> s^-1 m^-1 micron -> x1.E-6 -> s^-1
		"""

		photons_deltalambda = self.skyemission_photons * self.grating_trans \
					 * self.throughput	* self.area_tel * self.refl_tel
		return photons_deltalambda

	def compute_photonrate_spect(self):
		"""
		compute the number of target photons in a wavelength range determined by the selected grating/central wave
		INPUTS:
		"""
		"""
		   Units: [flux_band] is expected in J s^-1 m^-2 micron^-1
				  [h] = J s
				  [lambda_band,width_band] = micron
				  [telarea] = m^2
				  [throughput] = No units
					Results as ===> s^-1 m^-1 micron -> x1.E-6 -> s^-1
		"""
		photons_deltalambda = self.obj_photons_wave  
		photons_deltalambda *= self.atmostrans
		photons_deltalambda *= self.grating_trans
		photons_deltalambda *= self.throughput
		photons_deltalambda *= self.area_tel
		photons_deltalambda *= self.refl_tel
		return photons_deltalambda

	def signal_noise_texp_spect(self,Nexp_min, Nexp_max, Nexp_step,dit,aperture,nexp=1):
		"""
		Compute the S/N ratio as function of exposure time
		:param dit: detector integration time of individual exposure
		:param aperture:
		:param nexp=1: number of exposures
		:return:
		"""
		dit = dit * u.s
		nexp = np.arange(Nexp_min,Nexp_max,Nexp_step)
		texp = len(nexp) * dit
		#print ("nexp ",nexp[0],nexp[-1])
		print("aperture:",aperture)
		selected_index = 0			##### FIXME... Cual seleccionamos?
		phi_obj_apert = self.phi_obj_total * aperture['EE'][selected_index]
		phi_sky_apert = self.phi_sky_sangle * np.pi * aperture['Radius'][selected_index]**2
		print ("===========================")
		print (phi_obj_apert.unit) 
		print (phi_sky_apert.unit)
		print (self.detector['darkc']*aperture['Npix'][selected_index]).unit
		print (nexp*aperture['Npix'][selected_index]*self.detector['ron']**2).unit
		print (texp.unit)
		print ("===========================")
		noise = np.sqrt(texp.value*(phi_obj_apert.value+phi_sky_apert.value+(self.detector['darkc']*aperture['Npix'][selected_index]).value)+ (texp*aperture['Npix'][selected_index]*self.detector['ron']**2).value)
		signal = texp * phi_obj_apert
		print (signal.unit)
		snr = signal / noise
		return (texp,snr)

### end added by JAP

path_absolute = os.path.realpath(__file__).split("/")
path_absolute[len(path_absolute)-1] = ""
path_absolute="/".join(path_absolute)
if ( __name__ == "__main__"):

	conditions = Conditions(0.4)
	a = Calculator(0.0086,3.9e-6, 1, conditions, 1.64e-6, 1450)
	print('%lf' % a.test())
