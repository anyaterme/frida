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
	def __init__ (self, gs_mag=8, gs_dis=0.):
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
	def __init__ (self,mag_target,band,mag_system,sed,waveshift=('radial_velocity', 0.),\
	       wave_ini=8500*u.angstrom,wave_end=25000.*u.angstrom,dwave=2.*u.angstrom):
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
		self.sed_wave = sed_wave[mask]
		self.sed_flambda = sed_flambda[mask]


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
		wave_sameunit = wave.to(self.sed_wave.unit).value
		flambda = np.interp(wave_sameunit,self.sed_wave.value,self.sed_flambda.value)
		return flambda * flambda_unit

	def photons_wave(self,wave):
		"""
		Compute f_lambda for the selected SED and a given wave array
		:param wave:
		:return:
		"""
		energy_per_photon = const.h * const.c/ wave / u.photon
		phot_lambda = self.flambda_wave(wave)/(energy_per_photon) #const.h/const.c*wave
		return phot_lambda


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

	def __init__(self,target_info,filter,atmospheric_cond=None,scale=None,telescope_name='GTC',instrument_name='FRIDA'):
		"""
			Celestial source parameters in a dictionary
			target_info = {'Magnitude':19,'MagSystem':,'Band':'J','sed':sed} ;
			 sed is a tuple like ('PowerLaw',index_powerlaw), ('Blackbdoy',temperature_bb)
		"""
		sed = target_info.SED
		#temp_bb = sed[1]
		#input_band = target_info.Band   # definimos como banda de referencia
		#input_lambda= phot_zp[input_band][0]   # calculamos el flujo de referencia
		#mag_ref = float(target_info.Magnitude)
		#input_normflux=bbody(input_lambda,temp_bb)  # calculamos el flujo en la banda de referencia
		#flux_scale_old =scale_to_vega(input_band,mag_ref,input_normflux)
		flux_scale=target_info.flux_scale

		print ("flux scale",flux_scale,flux_scale)

		self.obs_filter_name=filter

		# telescope
		params_tel = Telescope(telescope_name)
		self.area_tel = params_tel.area
		self.refl_tel = params_tel.reflectivity
		print ('Telescope is ',params_tel.name)

		#### filter
		filter_info = Filter(filter,path_list=settings.INCLUDES,path_filters=settings.FILTERS)
		self.filter_wave =  filter_info.wave
		self.filter_trans = filter_info.transmission
		lambda_center = float(filter_info.lambda_center())
		filter_width = filter_info.width()
		print ("lambda_eff,bandwidth=",lambda_center,filter_width["width"],filter_width["cut_on"],filter_width["cut_off"])


		## instrument characteristics
		instrument=Instrument_static(scale,instid=instrument_name)
		static_response = instrument.compute_static_response(self.filter_wave)
		self.throughput = static_response["collimator"] * static_response['camera'] * static_response['qe']
		print ("Throughput ",self.throughput[1:5])
		self.detector = instrument.detector

		# Atmosphere
		atmosphere = Atmosphere()
		self.atmostrans = atmosphere.compute_skytrans(self.filter_wave)
		self.skyemission_photons = atmosphere.compute_skyrad(self.filter_wave)
		self.skymag = atmosphere.compute_skymag(self.obs_filter_name)

		# compute average atmospheric transmission
		atmospheric_trans_avg = atmosphere.compute_skytrans_avg(filter_info.wave[filter_width["index_above"]])
		print ("atmos_trans_avg=",atmospheric_trans_avg)

		collimator_transmission_avg = instrument.collimator_transmission.avg_response(filter_info.wave[filter_width["index_above"]])
		camera_transmission_avg = instrument.camera_transmission.avg_response(filter_info.wave[filter_width["index_above"]])
		detector_qe_avg = instrument.qe.avg_response(filter_info.wave[filter_width["index_above"]])
		filter_transmission_avg = filter_info.avg_transmission()
		print ("collimator_transmission_avg",collimator_transmission_avg)
		print ("camera_transmission_avg",camera_transmission_avg)
		print ("filter_transmission_avg",filter_transmission_avg)
		print ("detector_qe_avg",detector_qe_avg)
		print ("telescope reflectivity",self.refl_tel)
		self.effic_total = self.refl_tel * collimator_transmission_avg * \
				  filter_transmission_avg * camera_transmission_avg * detector_qe_avg
		print ("Instrument - Global efficiency",self.effic_total)



		## compute the target f-lambda within the wavelength range determined by the filter transmission
		obj_flambda_spec = target_info.flambda_wave(self.filter_wave)
		obj_photons_spec = target_info.photons_wave(self.filter_wave)
		ii = len(self.filter_wave)/2
		print (" wave,obj_flambda,photons["+str(ii)+"]=",self.filter_wave[ii],obj_flambda_spec[ii],\
			   obj_photons_spec[ii])


		## compute the flux from the target and the number of photons per unit time
		#self.obj_flambda = flux_scale * bbody(filter_info['lamb_eff'],temp_bb) # Units are
		self.obj_flambda_lambcenter = target_info.flambda_wave(lambda_center)
		print ("obj_flambda_lambeff:",self.obj_flambda_lambcenter)
		phi_obj_total2 = self.compute_photonrate_simple(self.obj_flambda_lambcenter,lambda_center,\
									   filter_width["width"],atmospheric_trans_avg*self.effic_total)

		self.flambda_wave = target_info.flambda_wave
		self.phi_obj_total = self.compute_photonrate_filter(target_info.flambda_wave(self.filter_wave))

		print("phi_obj_total [simple int-filter]",phi_obj_total2,self.phi_obj_total)


		## compute the photon rate from the sky per arcsec^2
		##sky_flambda = 10.**(phot_zp[self.obs_filter_name][2]-0.4*self.atmosphere['mag_sqarc'])
		sky_flambda = self.get_skyflux_band(lambda_center)
		phi_sky_sqarc_simple = self.compute_photonrate_simple(sky_flambda,lambda_center,\
									   filter_width["width"],self.effic_total)
		self.phi_sky_sqarc = self.compute_sky_photonrate_filter()

		print("phi_sky/pix/sec (simple  int-filter)",phi_sky_sqarc_simple*instrument.pixscale**2,\
			  self.phi_sky_sqarc*instrument.pixscale**2)


	def compute_photonrate_simple(self,flux_filter,lambda_filter,width_filter,throughput):
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
		phot_tot = flux_filter/const.h/const.c*lambda_filter*width_filter*self.area_tel*throughput  # units are photons/s
		return phot_tot

	def compute_photonrate_filter(self,flambda):
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
					Results as ===> s^-1 m^-1 micron -> x1.E-6 -> s^-1
		"""
		import scipy.integrate

		a = self.integrand_obj(flambda) # include convolution with atmospheric and filter transmission
		photons_through_filter = scipy.integrate.simps(a,self.filter_wave)
		return photons_through_filter*self.area_tel*self.refl_tel


	def integrand_obj(self, flux_obj):
		print ("Integrating Object Photon Flux")
		print (" Wave   Fl_obj  Atm_trans  Filt_trans  Phot_Ener   Throughput")
		ii = 150
		print (self.filter_wave[ii], flux_obj[ii], self.atmostrans[ii], self.filter_trans[ii],self.energy_photon(self.filter_wave[ii]), self.throughput[ii])
		#print ("%7.3F"% self.filter_wave[ii].value, "%9.2E"% flux_obj[ii].value, "%6.2F"% self.atmostrans[ii].value, "%6.2F"% self.filter_trans[ii].value,"%9.2E"% self.energy_photon(self.filter_wave[ii].value), "%6.2F"% self.throughput[ii].value)
		#print ("%7.3F"% self.filter_wave[ii].value, "%9.2E"% flux_obj[ii].value, "%6.2F"% self.atmostrans[ii].value, "%6.2F"% self.filter_trans[ii].value,"%9.2E"% self.energy_photon(self.filter_wave[ii].value), "%6.2F"% self.throughput[ii].value)

		return( flux_obj * self.atmostrans * self.filter_trans / self.energy_photon(self.filter_wave) * self.throughput)

	def energy_photon(self,wave):
		"""
		Compute the energy of a photon at a given wavelenght [in microns]
		:param wave:
		:return:
		"""
		return self.energy_photon_1mic/wave

	def compute_sky_photonrate_filter(self):
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
		import scipy.integrate

		a = self.integrand_sky() # include convolution with filter transmission
		photons_through_filter = scipy.integrate.simps(a,self.filter_wave)
		return photons_through_filter*self.area_tel*self.refl_tel

	def integrand_sky(self):
		return(self.skyemission_photons * self.filter_trans * self.throughput)

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
		#print ("nexp ",nexp[0],nexp[-1])
		phi_obj_apert = self.phi_obj_total * aperture['EE']
		phi_sky_apert = self.phi_sky_sqarc * pi * aperture['Radius']**2
		noise = np.sqrt(texp*(phi_obj_apert+phi_sky_apert+self.detector['darkc']*aperture['Npix'])+ \
			nexp*aperture['Npix']*self.detector['ron']**2)
		signal = texp * phi_obj_apert
		snr = signal / noise
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
		phi_sky_apert = self.phi_sky_sqarc * pi * aperture['Radius']**2
		noise_ndit = dit*(phi_obj_apert+phi_sky_apert+self.detector['darkc']*aperture['Npix'])+aperture['Npix']*self.detector['ron']**2

		ndit = round(required_sn**2 * noise_ndit / phi_obj_apert**2 / dit**2)

		return (ndit)


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
#		print "=============DEBUG=============", wave_zp, photzp_values, wave
#		phot_zp_interp = interpolate(np.array(wave_zp),np.array(photzp_values),wave)
#
#		skyflux = 10.**(phot_zp_interp-0.4*self.skymag)

		print("================= DEBUG ==================")
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

	def __init__(self,target_info,sky_conditions,spectrograph_setup,scale,Atmospheric_Cond=None,Telescope='GTC',Instrument='FRIDA'):
		"""
			Celestial source parameters in a dictionary
			target_info = {'Magnitude':19,'MagSystem':,'Band':'J','sed':sed} ;
			 sed is a tuple like ('PowerLaw',index_powerlaw), ('Blackbdoy',temperature_bb)
		"""
		# Setup the wavelength array from spectrograph_setup
		#### grating
		grating_info = Grating(spectrograph_setup['Grating_name'],cent_wave=spectrograph_setup['Wave_center'])
		self.grating_wave =  grating_info.wave
		self.grating_trans = grating_info.transmission
		self.grating_width = grating_info.cutoff-grating_info.cuton
		print ("Disp,lambda_center,bandwidth=",grating_info.dlambda,grating_info.lcenter,self.grating_width)
		self.wave_array = grating_info.wave_array()
		self.grating_effic = grating_info.interp_spect("EfficGrating")

		sed = target_info.SED
		#temp_bb = sed[1]
		#input_band = target_info.Band   # definimos como banda de referencia
		#input_lambda= phot_zp[input_band][0]   # calculamos el flujo de referencia
		#mag_ref = float(target_info.Magnitude)
		#input_normflux=bbody(input_lambda,temp_bb)  # calculamos el flujo en la banda de referencia
		#flux_scale_old =scale_to_vega(input_band,mag_ref,input_normflux)
		flux_scale=target_info.flux_scale
		print ("flux scale=",flux_scale)
		photons_obj=target_info.photons_wave(self.wave_array)

		# determine flux density for the

		# telescope
		params_tel = Telescope(Telescope)
		self.area_tel = params_tel.area
		self.refl_tel = params_tel.reflectivity
		print ('Telescope is ',params_tel.name)


		## instrument characteristics  instid by default is FRIDA, other options are NACO
		instrument=Instrument_static(scale,instid=Instrument)
		static_response = instrument.compute_static_response(self.wave_array)
		self.throughput = static_response["collimator"] * static_response['camera'] * static_response['qe']
		print ("Throughput ",self.throughput[1:5])
		self.detector = instrument.detector

		# Atmosphere
		atmosphere = Atmosphere()
		## compute atmospheric transmission and emission, convolved with the expected resolution
		self.atmostrans = atmosphere.compute_skytrans(self.wave_array,kernel='gauss')
		self.skyemission_photons = atmosphere.compute_skyrad(self.wave_array,kernel='gauss')


		## compute the target f-lambda within the wavelength range determined by the filter transmission
		self.obj_flambda_wave = target_info.flambda_wave(self.wave_array)
		self.obj_photons_wave = target_info.photons_wave(self.wave_array)
		ii = len(self.wave_array)/2
		print (" wave,obj_flambda,photons["+str(ii)+"]=",self.wave_array[ii],obj_flambda_wave[ii],\
			   obj_photons_wave[ii])

		self.flambda_wave = target_info.flambda_wave
		# compute the effective photon rate from target and sky as measured by the detector
		self.phi_obj_total = self.compute_photonrate_spect(self)
		self.phi_sky_sangle = self.compute_sky_photonrate_spect(self)


		## compute the photon rate from the sky per arcsec^2
		##sky_flambda = 10.**(phot_zp[self.obs_filter_name][2]-0.4*self.atmosphere['mag_sqarc'])
		sky_flambda = self.get_skyflux_band(lambda_center)
		self.phi_sky_sqarc = self.compute_photonrate_simple(sky_flambda,lambda_center,\
									   filter_width["width"],self.effic_total)

		print("phi_sky/pix/sec (simple  int-filter)",self.phi_sky_sqarc*instrument.pixscale**2,\
			  self.compute_sky_photonrate_filter()*instrument.pixscale**2)

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

		photons_deltalambda = self.obj_photons_wave * self.atmostrans * self.grating_trans \
			         * self.throughput * self.area_tel*self.refl_tel
		return photons_deltalambda

	def signal_noise_texp_spect(self,dit,aperture,nexp=1):
		"""
		Compute the S/N ratio as function of exposure time
		:param dit: detector integration time of individual exposure
		:param aperture:
		:param nexp=1: number of exposures
		:return:
		"""
		nexp = np.arange(Nexp_min,Nexp_max,Nexp_step)
		texp = nexp * dit
		#print ("nexp ",nexp[0],nexp[-1])
		phi_obj_apert = self.phi_obj_total * aperture['EE']
		phi_sky_apert = self.phi_sky_sangle * pi * aperture['Radius']**2
		noise = np.sqrt(texp*(phi_obj_apert+phi_sky_apert+self.detector['darkc']*aperture['Npix'])+ \
			nexp*aperture['Npix']*self.detector['ron']**2)
		signal = texp * phi_obj_apert
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
