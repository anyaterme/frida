
from django.conf import settings
import astropy.units as u
import astropy.constants as const
import scipy.integrate
import numpy as np
import csv
import os
from frida.aux_functions import convolve2pulse,convolve2gauss
from frida.aux_functions import interpol2newwave,interpolate
from frida.aux_functions import read_grating_files

class Telescope:

	def __init__ (self, telid="GTC"):

		telarea = None
		refl_total = None
		teldiam1 = None

		if (telid == 'GTC'):
			name = '10m GTC'
			teldiam1=10.4 * u.m # in meter
			teldiam2=1. * u.m
			#telarea=3.14*(teldiam1*teldiam1-teldiam2*teldiam2)/4.
			telarea=73.4  * u.m**2   ## in m^2 (provided by A. Cabrera-GTC)
			refl_primary=0.9
			refl_second=0.9
			refl_tertia=0.9
			refl_total=refl_primary*refl_second*refl_tertia
		elif (telid == 'VLT'):
			name = '8m VLT'
			teldiam1=8.2 * u.m # in meter
			teldiam2=1.* u.m
			telarea=3.14*(teldiam1*teldiam1-teldiam2*teldiam2)/4.
			refl_primary=0.9
			refl_second=0.9
			refl_total=refl_primary*refl_second
		elif (telid == 'WHT'):
			name = '4.2m WHT'
			teldiam1=4.2 * u.m # in meter
			teldiam2=1.* u.m
			telarea=3.14*(teldiam1*teldiam1-teldiam2*teldiam2)/4.
			refl_primary=0.9
			refl_second=0.9
			refl_total=refl_primary*refl_second

		self.name = name
		self.area = telarea
		self.reflectivity = refl_total
		self.aperture = teldiam1

class Spec_response:
	def __init__(self,file,path=settings.INCLUDES):
		#fp = np.loadtxt(os.path.join(path,file))
		#self.wave = fp[:,0] * u.micron
		#self.response = fp[:,1]
		fp = open(os.path.join(path, file))
		sprf_file = csv.DictReader(filter(lambda row: row[0] != '#', fp),delimiter=' ',skipinitialspace=True)
		sprf_fn0 = sprf_file.fieldnames[0]
		sprf_fn1 = sprf_file.fieldnames[1]
		sprf_wave = []
		sprf_value = []
		for row in sprf_file:
			sprf_wave.append(float(row[sprf_fn0]))
			sprf_value.append(float(row[sprf_fn1]))
		sprf_wave = np.array(sprf_wave)
		sprf_value = np.array(sprf_value)
		fp.close()
		# give units
		if '[nm]' in sprf_fn0:
			sprf_wave = sprf_wave * u.nm
		if '[angstrom]' in sprf_fn0:
			sprf_wave = sprf_wave * u.AA
		if '[micron]' in sprf_fn0:
			sprf_wave = sprf_wave* u.micron
		## if value is given in % convert to fraction
		if '[%]' in sprf_fn1:
			sprf_value = sprf_value / 100.
		self.wave = sprf_wave
		self.response = sprf_value

	def interpol2wave(self,wave):
		#return np.interp(wave,self.wave,self.response)
		return np.interp(wave,self.wave,self.value)

	def avg_response(self,wave_range):
		# first interpolate
		#response_wave =  interpolate(self.wave,self.response,wave_range)
		response_wave =  interpolate(self.wave,self.value,wave_range)
		avg_response = scipy.integrate.simps(response_wave,wave_range)/(wave_range[-1]-wave_range[0])
		return avg_response

class Instrument_static:
	"""
	This module will load all parameters related to the instrument, e.g. the transmission files, the detector characteristics
	"""
	def __init__ (self,scale='fine_scale',instid='FRIDA',path=settings.INCLUDES):
		self.qe = Spec_response("detector_qe.dat",path=path)
		self.qe.response = self.qe.response * u.electron/u.photon

		if (scale == 'fine_scale'):
			self.pixscale = 0.010 * u.arcsec
			if (instid == "NACO"):
				self.pixscale = 0.013 * u.arcsec
			self.camera_trans = Spec_response("finecamera_transmission.dat",path=path)
			print ("Selecting fine scale")
			print ("camera _trans=",self.camera_trans.response[0:3])
		elif (scale == 'medium_scale'):
			self.pixscale = 0.02 * u.arcsec
			self.camera_trans = Spec_response("mediumcamera_transmission.dat",path=path)
		elif (scale == 'coarse_scale'):
			self.pixscale = 0.04 * u.arcsec
			self.camera_trans = Spec_response("coarsecamera_transmission.dat",path=path)

		self.collimator_trans = Spec_response("collimator_transmission.dat",path=path)

		## FRIDA Detector parameters
		read_out_noise = 15. * u.electron   ## in e-
		gain = 3.5  # e-/ADU
		welldepth = 2.e5 * u.electron ## in e-
		dark_current = 0.1 * u.electron / u.second ## dark current [e-/sec]

		self.detector = {"ron":read_out_noise,"gain":gain,"welldepth":welldepth,"darkc":dark_current}

	def compute_static_response(self,nwave):
		"""
		Interpolate all values of qe, collimator and camera transmission
		:param wave:
		:return:
		"""
		intp_qe = interpol2newwave(self.qe.response,self.qe.wave,nwave)
		intp_camera = interpol2newwave(self.camera_trans.response,self.camera_trans.wave,nwave)
		intp_collim = interpol2newwave(self.collimator_trans.response,self.collimator_trans.wave,nwave)
		return {"qe":intp_qe,"collimator":intp_collim,"camera":intp_camera}

def param_gratings():
	#
	#       Gets the characteristics of filter: cut and returns it as a
	#       speccurve object
	#
	import csv
	fp = open(os.path.join(settings.INCLUDES, 'gratings.dat'))
	list_filters= csv.DictReader(filter(lambda row: row[0] != '#', fp))
	for each_filter in list_filters:
		if each_filter["Name"] == filter_name :
			myfilter = each_filter
			print("Filter Code: ", myfilter["Code"])
			print("Filter_Transmission File: ",myfilter["Transmission"])
			print("Lambda center",myfilter["lambda_center"])
			break
	fp.close()

	lambda_eff = float(myfilter["lambda_center"])
	cut_on = float(myfilter["cut-on"])
	cut_off = float(myfilter["cut-off"])

	bandwidth=(cut_off-cut_on)
	energy_phot = const.h*const.c/(lambda_eff*1.e-6) ## in SI units
	avg_trans = 0.8

	return {'name':myfilter["Name"],'lamb_eff':lambda_eff,'bandwidth':bandwidth,'avg-transmission':avg_trans,
			'energy_phot':energy_phot}


def param_filter(filter_name='H'):
	#
	#       Gets the characteristics of filter: cut and returns it as a
	#       speccurve object
	#
	import csv
	fp = open(os.path.join(settings.INCLUDES, 'filters.dat'))
	list_filters= csv.DictReader(filter(lambda row: row[0] != '#', fp))
	for each_filter in list_filters:
		if each_filter["Name"] == filter_name :
			myfilter = each_filter
			print("Filter Code: ", myfilter["Code"])
			print("Filter_Transmission File: ",myfilter["Transmission"])
			print("Lambda center",myfilter["lambda_center"])
			break
	fp.close()

	lambda_eff = float(myfilter["lambda_center"])
	cut_on = float(myfilter["cut-on"])
	cut_off = float(myfilter["cut-off"])

	bandwidth=(cut_off-cut_on)
	energy_phot = const.h*const.c/(lambda_eff*1.e-6) ## in SI units
	avg_trans = 0.8

	return {'name':myfilter["Name"],'lamb_eff':lambda_eff,'bandwidth':bandwidth,'avg-transmission':avg_trans,
			'energy_phot':energy_phot}


class Grating:

	def __init__ (self,grating_name,cent_wave=None,path_list=settings.INCLUDES,path_gratings=settings.GRATINGS,path_skycalc=settings.SKYCALC_GRATINGS):
		import csv
		fp = open(os.path.join(path_list, 'gratings.dat'))
		list_gratings= csv.DictReader(filter(lambda row: row[0] != '#', fp))
		for each_grating in list_gratings:
			print("Found grating  ",each_grating["Name"])
			if each_grating["Name"] == grating_name :
				mygrating = each_grating
				print("Grating Name: ", mygrating["Name"])
				print("Filter_Transmission File: ",mygrating["Efficiency"])
				print("Central Wavelength",mygrating["Central_Wave"])
				break
		fp.close()

		fp = open(os.path.join(path_list, 'gratings2skycalc.dat'))
		list_gratings= csv.DictReader(filter(lambda row: row[0] != '#', fp))
		for each_grating in list_gratings:
			print("Found grating  ",each_grating["Name"])
			if each_grating["Name"] == grating_name :
				mygratmosph = each_grating
				print("Grating Name: ", mygratmosph["Name"])
				print("Sky Radiance File: ",mygratmosph["Sky_radiance"])
				print("Sky Transmission: ",mygratmosph["Atmos_trans_WV2p5"])
				break
		fp.close()

		#ftrans = np.loadtxt(os.path.join(path_gratings, mygrating["Efficiency"]))
		## read efficiency file
		greffic=read_grating_files(mygrating["Efficiency"],path_gratings=path_gratings)
		## create an array using dispersion in the whole range available in the transmission curve
		self.delt_wave = float(mygrating["Dispersion"]) * u.AA
		if cent_wave is not None:
			self.cent_wave = cent_wave
		else:
			self.cent_wave = float(mygrating["Central_Wave"]) * u.AA
		self.cuton=mygrating["Wave_ini"]
		self.cutoff=mygrating["Wave_end"]
		self.gr_wave = greffic['wave']
		self.gr_effic = greffic['effic']
		print ("Wave efficiency:",self.gr_wave[0:4],self.gr_effic[0:4])

		# read sky emission and atmospheric absorption
		fp = open(os.path.join(path_skycalc, mygratmosph["Sky_radiance"]))
		skyrad_file = csv.DictReader(filter(lambda row: row[0] != '#', fp),delimiter=' ')
		skyrad_fn0 = skyrad_file.fieldnames[0]
		skyrad_fn1 = skyrad_file.fieldnames[1]
		skyrad_wave = []
		skyrad_value = []
		for row in skyrad_file:
			skyrad_wave.append(float(row[skyrad_fn0]))
			skyrad_value.append(float(row[skyrad_fn1]))
		skyrad_wave = np.array(skyrad_wave)
		skyrad_value = np.array(skyrad_value)
		# give units
		if '[nm]' in skyrad_fn0:
			skyrad_wave = skyrad_wave * u.nm
		if '[angstrom]' in skyrad_fn0:
			skyrad_wave = skyrad_wave * u.AA
		if '[micron]' in skyrad_fn0:
			skyrad_wave = skyrad_wave * u.micron
		## if value is given in % convert to fraction
		skyrad_value = skyrad_value * u.photon / u.s / u.m**2 / u.micron / u.arcsec**2
		self.skyrad = {'wave':skyrad_wave,'value':skyrad_value}
		fp.close()

# read sky emission and atmospheric absorption
		fp = open(os.path.join(path_skycalc, mygratmosph["Atmos_trans_WV2p5"]))
		atrans_file = csv.DictReader(filter(lambda row: row[0] != '#', fp),delimiter=' ')
		atrans_fn0 = atrans_file.fieldnames[0]
		atrans_fn1 = atrans_file.fieldnames[1]
		atrans_wave = []
		atrans_value = []
		for row in atrans_file:
			atrans_wave.append(float(row[atrans_fn0]))
			atrans_value.append(float(row[atrans_fn1]))
		atrans_wave = np.array(atrans_wave)
		atrans_value = np.array(atrans_value)
		# give units
		if '[nm]' in atrans_fn0:
			atrans_wave = atrans_wave * u.nm
		if '[angstrom]' in atrans_fn0:
			atrans_wave = atrans_wave * u.AA
		if '[micron]' in atrans_fn0:
			atrans_wave = atrans_wave * u.micron
		## if value is given in % convert to fraction
		self.atrans = {'wave':atrans_wave,'value':atrans_value}
		fp.close()

	def wave_array(self,npix=2048):
		"""
		This function generates an array with the wavelength reference to compute S/N in spectroscopy mode
		:param self.cent_wave:
		:param npix:
		:return:
		"""
		wave = (np.arange(0,npix,1)-npix/2)*self.delt_wave + self.cent_wave
		#wave = (np.arange(0,npix,1)-npix/2)*self.delt_wave
		return wave

	def interp_spect(self,curve_type,wave_new):
		"""
		This function generates an array with the wavele to compute S/N in spectroscopy mode
		:param wave_center:
		:param npix:
		:return:
		"""
		#wave_new=self.wave_array()
		#wave_int=self.gr_wave
		#wave_new = [1.,2.]
		#wave_int = [1.,2,.,4.]
		#interp_spect = interpolate(self.gr_wave, self.gr_effic, wave_array(self,npix=npix))
		###valint = {
		###	'GratingEffic': lambda x: np.interp(x,self.gr_wave.to('micron'),self.gr_effic),
		###	'AtmosTrans': lambda x: np.interp(x, self.atrans['wave'].to('micron'), self.atrans['value']),
		###	'SkyRad': lambda x: np.interp(x, self.skyrad['wave'].to('micron'), self.skyrad['value'])
		###} [curve_type](wave_new.to('micron'))
		valp = {'GratingEffic':self.gr_effic,'AtmosTrans':self.atrans['value'],'SkyRad':self.skyrad['value']}[curve_type]
		wavep = {'GratingEffic':self.gr_wave,'AtmosTrans':self.atrans['wave'],'SkyRad':self.skyrad['wave']}[curve_type]
		valunit = {'GratingEffic':1,'AtmosTrans':1,'SkyRad':self.skyrad['value'].unit}[curve_type]
		valint = np.interp(wave_new.to("micron"),wavep.to("micron"),valp)
		valint = valint * valunit
		return valint

class Filter:
	"""
	This class includes any processing to be done with filters in imaging mode
       wave - array containing the wavelength range around the 
              transmision range
       transmission - array containing the transmission
	"""
	def __init__ (self,filter_code,path_list=None,path_filters=None):
		import csv
		fp = open(os.path.join(path_list, 'filters.dat'))
		list_filters= csv.DictReader(filter(lambda row: row[0] != '#', fp))
		myfilter = None
		print ("Filter Name:", filter_code)
		for each_filter in list_filters:
			if myfilter is None:
				if os.path.isfile(os.path.join(path_filters, each_filter["Transmission"])):
					myfilter = each_filter
			if each_filter["Code"] == filter_code :
				myfilter = each_filter
				print("Filter Code: ", myfilter["Code"])
				print("Filter_Transmission File: ",myfilter["Transmission"])
				print("Lambda center",myfilter["lambda_center"])
				break
		fp.close()

		if (myfilter is None):
			print("Filter is not selected")
			##FIXME Add Default Values
		else:
			ftrans = np.loadtxt(os.path.join(path_filters, myfilter["Transmission"]))
			#ftrans = np.sort(ftrans,axis=0)
			ftrans=np.asarray(sorted(ftrans, key=lambda row: row[0]))
			self.wave = ftrans[:,0] * u.micron
			self.info = myfilter
			self.transmission = ftrans[:,1] * 0.01
			#self.wave_limits = (min(self.wave),max(self.wave))
			self.wave_limits = (float(myfilter["cut-on"]) *0.95 * u.micron, float(myfilter["cut-off"]) * 1.05 * u.micron) 
			self.wave_median = self.wave_limits[0] + np.abs(self.wave_limits[1]-self.wave_limits[0])/2.
			print ("Filter transmision:",self.wave[0:4],self.transmission[0:4])
			print ("Filter min, max:",self.wave_limits[0],self.wave_limits[1])

	def avg_transmission(self,threshold_fraction=2.):
		'''
          Compute the average transmission value over a threshold
		'''         
		wave_value = self.wave.value
		trans = self.transmission
		# compute transmision above threshold
		threshold_trans = trans.max()/threshold_fraction
		mask_above = [trans  >= threshold_trans]
		avg_trans = scipy.integrate.simps(trans[mask_above],wave_value[mask_above])\
             /(wave_value[mask_above][-1]-wave_value[mask_above][0])
		return avg_trans

	def lambda_center(self):
		'''
          Compute the effective wavelength for a given filter transmission curve
		'''         
		wave_unit = self.wave.unit
		wave_value = self.wave.value
		lambda_eff = scipy.integrate.simps(self.transmission*wave_value,wave_value)/\
            scipy.integrate.simps(self.transmission,wave_value)
		return lambda_eff * wave_unit

	def width(self,threshold_fraction=2.):
		trans = self.transmission
		threshold_trans = trans.max()/threshold_fraction
		index_above = np.where(trans >= threshold_trans)
		wave_above = self.wave[index_above]
		cut_on = wave_above[0]
		cut_off = wave_above[-1]
		print ("cut-on,off=",cut_on,cut_off)
		width = cut_off - cut_on
		return {"width":width,"cut_on":cut_on,"cut_off":cut_off,"index_above":index_above}



class Atmosphere:
	"""
	This class includes any processing dealing with Earth atmosphere transmission and emission
	"""

	def __init__ (self, transmission="skytransmission_mean.dat",radiance="skyradiance_mean.dat",path=settings.INCLUDES):
		f = np.loadtxt(os.path.join(path,transmission))
		## Wavelenght units are given in nm
		self.atmtrans_wave = f[:,0] * u.nm
		self.atmtrans_perone = f[:,1]
		self.atmtrans_delta = np.abs(self.atmtrans_wave[1]-self.atmtrans_wave[0])
		print("atmtrans_delta:",self.atmtrans_delta)

		f = np.loadtxt(os.path.join(path,radiance))
		## Wavelenght units are nm, must be converted to microns
		self.skyrad_wave = f[:,0] * u.nm
		self.skyrad_photons = f[:,1] * u.photon / u.s/(u.m)**2/u.micron/(u.arcsec)**2 
		self.skyrad_delta = np.abs(self.skyrad_wave[1]-self.skyrad_wave[0])
		print("skyrad_delta:",self.skyrad_delta)

	def compute_skytrans(self,wave,kernel='gauss'):
		# first check spacing of wave
		delta = abs(wave[1]-wave[0])
		#delta = delta.value
		print ("  @compute_skytrans - delta, atmtrans_delta",delta,self.atmtrans_delta)
		if (self.atmtrans_delta < delta):
			if (kernel == 'pulse'):
				print("  @compute_skytrans - Using convolve pulse")
				atmtrans_conv = convolve2pulse(self.atmtrans_wave,self.atmtrans_perone,2*delta)
			elif (kernel == 'gauss'):
				print("  @compute_skytrans - Using convolve gauss")
				atmtrans_conv = convolve2gauss(self.atmtrans_wave,self.atmtrans_perone,2*delta)
			atmtrans_interp = interpolate(self.atmtrans_wave,atmtrans_conv,wave)
		else:
			atmtrans_interp = interpolate(self.atmtrans_wave,self.atmtrans_perone,wave)
		print(atmtrans_interp)
		return atmtrans_interp

	def compute_skyrad(self,wave,kernel='pulse'):
		# first check spacing of wave
		delta = abs(wave[1]-wave[0])
		print ("  @compute_skyrad - delta, atmtrans_delta",delta,self.skyrad_delta)
		if (self.skyrad_delta < delta):
			if (kernel == 'pulse'):
				skyrad_conv = convolve2pulse(self.skyrad_wave,self.skyrad_photons,2*delta)
			elif (kernel == 'gauss'):
				skyrad_conv = convolve2gauss(self.skyrad_wave,self.skyrad_photons,2*delta)
			skyrad_interp = interpolate(self.skyrad_wave,skyrad_conv,wave)
		else:
			skyrad_interp = interpolate(self.skyrad_wave,self.skyrad_photons,wave)
		return skyrad_interp

	def compute_skytrans_avg(self,wave):
		# compute average sky transmision in a wavelength range
		sky_trans = self.compute_skytrans(wave)
		skytrans_avg = scipy.integrate.simps(sky_trans,wave)/(wave[-1]-wave[0])
		return skytrans_avg

	def compute_skymag(self,filter='H'):
		skymag = 13.
		if (filter == 'H'):
			skymag=14.4
		elif (filter == 'J'):
			skymag=17.
		elif (filter == 'Ks'):
			skymag=13.

		return skymag


def flux_vega(filter='H'):

	if (filter == 'H'):
		name='H'
		f0=1.14e-9
	elif (filter == 'J'):
		name='J'
		f0=3.15e-9
	elif (filter == 'Ks'):
		name='Ks'
		f0=3.96e-10
	elif (filter == 'Jc'):
		name='Jc'
		f0=1.14e-9

	units='W m-2 mic-1'
	#
	return {'name':name,'flux0':f0,'units':units}


path_absolute = os.path.realpath(__file__).split("/")
path_absolute[len(path_absolute)-1] = ""
path_absolute="/".join(path_absolute)
