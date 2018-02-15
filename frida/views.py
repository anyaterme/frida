"""
A view function, or view for short, is simply a Python function that takes a Web request and returns a Web response.
This response can be the HTML contents of a Web page, or a redirect, or a 404 error, or an XML document, or
an image . . . or anything, really.
This module reads information from web  (params1,2 and 3) and redirects to

"""
from django.shortcuts import render
from django.conf import settings
from django.core.files import File
from django.http import HttpResponse
import os
import datetime
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#from math import sqrt, pi, log10

#import frida.calculator
#import frida.gtcao
#import frida.set_config
from frida.set_config import *
from frida.compute_flux import define_photzp
from frida.calculator import SkyConditions, GuideStar, TargetInfo, Calculator_Image, Calculator_Spect
from frida.gtcao import GTC_AO
#from django.http import JsonResponse

# Create your views here.
def index(request):
	import csv

	#fp.close()
	#list_filters = {}

	fp1 = open(os.path.join(settings.INCLUDES, 'filters.dat'))
	list_filters = csv.DictReader(filter(lambda row: row[0] != '#', fp1))
	list_filters_dict = []
	for each_filter in list_filters:
		my_filter = {}
		my_filter["Code"] = each_filter["Code"]
		my_filter["Name"] = each_filter["Name"]
		list_filters_dict.append(my_filter)

	#fp1.close()

	fp2 = open(os.path.join(settings.INCLUDES, 'gratings.dat'))
	list_gratings= csv.DictReader(filter(lambda row: row[0] != '#', fp2))
	list_gratings_dict = []
	for each_grating in list_gratings:
		my_grating = {}
		my_grating["Description"] = each_grating["Description"]
		my_grating["Name"] = each_grating["Name"]
		my_grating["Central_Wave"] = each_grating["Central_Wave"]
		my_grating["Efficiency"] = each_grating["Efficiency"]
		list_gratings_dict.append(my_grating)
	#fp2.close()

	context = {'list_filters':list_filters_dict, 'list_gratings':list_gratings_dict}
	return render(request, 'index.html', context)

def calculate_draw(request):
	context = {}
	return render(request, 'result.html', context)


def calculate(request):
	"""
	select according to observing mode:
		-imaging -> calculate_ima
		-IFS -> calculate_ifs
	:param request:
	:return:
	"""
	if (request.POST.get('observation_mode') == 'IFS'):
		return (calculate_ifs(request))
	else:
		return (calculate_ima(request))

def get_SkyConditions(request):
	"""
	Obtain values about Sky conditions
	:param request: input from html
	:return: dict sky_conditions{'seeing','lambda_seeing','airmass','pwv'}
	"""
	# cond = Conditions()
	sky_conditions = {'seeing': 0.9, 'lambda_seeing': 0.5 * u.micron, 'airmass': 1.2, 'pwv': 2.5}
	sky_conditions['seeing'] = float(request.POST.get("seeing", sky_conditions['seeing']))
	sky_conditions['airmass'] = float(request.POST.get("airmass", sky_conditions['airmass']))
	# (JAP) FALTA INCLUIR VAPOR DE AGUA
	return sky_conditions

def get_TargetInfo(request):
	"""
	Obtain requirements about target from html templates. It returns an object of class TargetInfo
	:param request: input from html
	:return: object TargetInfo
	"""
	debug_values = {}
	# Request parameters relative to astronomical target
	mag_target = float(request.POST.get('spatial_integrated_brightness'))
	band = request.POST.get("band_flux")
	mag_system = request.POST.get('units_sib')
	energy_type = request.POST.get('spectral_type')
	debug_values["Brightness (Astronomical Source Definition)"] = mag_target
	debug_values["Filter_Name (Astronomical Source Definition)"] = band
	debug_values["Energy Distribution (Astronomical Source Definition)"] = energy_type

	## Select morphology of source
	source_morpho = request.POST.get('source_type')
	if (source_morpho == 'extended'):
		label_source_morpho = 'Extended source'
	elif (source_morpho == 'point'):
		label_source_morpho = 'Point source'

	# create a tuple with the information relative to the Spectral Energy Distribution
	if (energy_type == 'black_body'):
		temperature= float(request.POST.get('temperature'))
		sed = (energy_type,temperature)
		label_energy_type = 'Black Body, T=%s K' % temperature
	elif (energy_type == 'power_law'):
		pl_index = float(request.POST.get('pl_index'))
		sed = (energy_type,pl_index)
		label_energy_type = 'Power Law, lambda^%s' % pl_index

	# creates an object of type TargetInfo, it provides method to compute scaled f-lambda at any wavelength
	target_info = TargetInfo(mag_target,band,mag_system,sed)
	return target_info

def calculate_ima(request):
	"""
    Calculations for imaging mode. It reads parameters from the target, the sky conditions (), the instrument setup
    (pixel scale,  filter).
    Depending on the choice it computes S/N ratio for a range of exposure times or determine the exposure time to reach
    a desired S/N ratio.

	:param request:
	:return:
    """

	debug_values = {}

	sky_conditions = get_SkyConditions(request)
	"""
	#cond = Conditions()
	sky_conditions={'seeing':0.9,'lambda_seeing':0.5,'airmass':1.2,'pwv':2.5}
	try:
		sky_conditions['seeing'] = float(request.POST.get("seeing"))
		sky_conditions['airmass'] = float(request.POST.get("airmass"))
	except:
		pass
    """

	selected_scale = request.POST.get('scale')
	frida_setup = Instrument_static(selected_scale,path=settings.INCLUDES)

	selected_filter = request.POST.get('filter')
	msg_err = []
	try:
		obs_filter = Filter(selected_filter,path_list=settings.INCLUDES,path_filters=settings.FILTERS)
	except Exception as e:
		msg_err.append("%s" % e)
		return render(request, "error.html", {'msg_err':msg_err})
		
	lambda_eff=obs_filter.lambda_center()


	# creates an object of type TargetInfo, it provides method to compute scaled f-lambda at any wavelength
	#target_info = TargetInfo(mag_target,band,mag_system,sed)
	target_info = get_TargetInfo(request)
	"""  OLD WAY
	# Request parameters relative to astronomical target
	mag_target = float(request.POST.get('spatial_integrated_brightness'))
	band = request.POST.get("band_flux")
	mag_system = request.POST.get('units_sib')
	energy_type = request.POST.get('spectral_type')
	debug_values["Brightness (Astronomical Source Definition)"] = mag_target
	debug_values["Filter_Name (Astronomical Source Definition)"] = band
	debug_values["Energy Distribution (Astronomical Source Definition)"] = energy_type

	## Select geometry of source
	source_type = request.POST.get('source_type')
	if (source_type == 'extended'):
		label_source_type = 'Extended source'
	elif (source_type == 'point'):
		label_source_type = 'Point source'

	# create a tuple with the information relative to the Spectral Energy Distribution
	if (energy_type == 'black_body'):
		temperature= float(request.POST.get('temperature'))
		sed = (energy_type,temperature)
		label_energy_type = 'Black Body, T=%s K' % temperature
	elif (energy_type == 'power_law'):
		pl_index = float(request.POST.get('pl_index'))
		sed = (energy_type,pl_index)
		label_energy_type = 'Power Law, lambda^%s' % pl_index
    """

	telescope = "VLT"
	telescope_params = Telescope(telescope)

	# Creates an object of type Calculator_Image
	a=Calculator_Image(target_info,selected_filter,sky_conditions,selected_scale,telescope_name=telescope)


    ## call GTC_AO to compute Strehl ratio and Encircled Energy
	guide_star = GuideStar(float(request.POST.get("gs_magnitude")),float(request.POST.get("gs_separation")))

	aocor = GTC_AO(sky_conditions,guide_star,telescope_params.aperture)
	strehl= aocor.compute_strehl(lambda_eff,sky_conditions['airmass'])
	psf = aocor.compute_psf(lambda_eff,strehl['StrehlR'])
	fcore = 1.5 # radius of aperture as a factor of FWHM of the PSF core
	pixscale = frida_setup.pixscale # in arcseconds
	ee_aperture=aocor.compute_ee(psf,pixscale,fcore=fcore)

	calc_method = request.POST.get('type_results')

	dit = float(request.POST.get('DIT_exp'))
	Nexp = int(request.POST.get('N_exp'))
	Nexp_min = 1
	Nexp_max = 100.
	Nexp_step = (Nexp_max-Nexp_min)/30
	(texp_seq,snr_seq) = a.signal_noise_texp_img(Nexp_min,Nexp_max,Nexp_step,dit,ee_aperture)

	signal_noise= float(request.POST.get('signal_noise'))
	if (calc_method=="SN_ratio"):
		required_sn = signal_noise
	else:
		required_sn = snr_seq[-1]
	ndit = a.texp_signal_noise_img(required_sn,dit,ee_aperture)
	print (" Required S/N ",required_sn)
	print ("Time required to reach S/N ",ndit,ndit*dit)

	################################
	# Graphical output
	# exposure time vs SNR
	len_texp_seq = len(texp_seq)
	m4 = np.zeros((len_texp_seq,2))
	print ("texp, snr ",texp_seq[len_texp_seq-1],snr_seq[len_texp_seq-1])
	m4[:,0] = texp_seq
	m4[:,1] = snr_seq
	m4 = str(m4.tolist())

	# wavelength vs SED above atmosphere
	nele = len(obs_filter.wave)
	m2 = np.zeros((nele,2))
	m2[:,0] = obs_filter.wave
	m2[:,1] = target_info.flambda_wave(obs_filter.wave)/1.e-16
	m2 = str(m2.tolist())


	debug_values["FWHM_core"] = psf['FWHM_core']

	context = {
					'graph_title':'Image Mode',
					'sed_flambda':m2,
					'energy_type': label_energy_type,
					'source_type': label_source_type,
					'sky_conditions' : sky_conditions,
					'target_info': target_info,
					'guide_star': guide_star,
					'frida_setup': frida_setup,
					'filter': obs_filter,
					'signal_noise_req': required_sn,
					'total_exposure_time': dit*ndit,
					'signal_noise':m4,
					'lambda_center': obs_filter.lambda_center(),
					'dit' : dit,
					'ndit' : ndit,
					'Object_magnitude': mag_target,
					'Input_Band': band,
					'strehl_ratio':strehl['StrehlR'],
					'encircled_energy': ee['EE'],
					'Aperture_radius': aperture['radius'],
					'Pixscale': pixscale,
					'AreaNpix': aperture['Npix'],
					'flux_obj': a.obj_flambda_lambcenter,
					'detected_photons_from_source': str("%9.2E"% a.phi_obj_total),
					'detected_photons_from_sky_sqarcsec' : str("%9.2E"% a.phi_sky_sqarc),
					'atmospheric_transmission': a.atmostrans[10],
					'efficiency': a.effic_total,
					'debug_values': debug_values,
					}
	return render(request, 'result.html', context)


def calculate_ifs(request,telescope=settings.TELESCOPE):
	debug_values = ["test ==> test_content"]

	sky_conditions = get_SkyConditions(request)

	telescope_params = Telescope(telescope)

	selected_scale = request.POST.get('scale')
	frida_setup = Instrument_static(scale=selected_scale,instid=settings.INSTRUMENT,path=settings.INCLUDES)

	selected_grating = request.POST.get('grating')
	central_wave_grating = request.POST.get('grating_value',None)
	spectrograph_setup = {'Grating':selected_grating,'Central_wave':float(central_wave_grating)*u.AA}

	# creates an object of type TargetInfo, it provides method to compute scaled f-lambda at any wavelength
	#target_info = TargetInfo(mag_target,band,mag_system,sed)
	target_info = get_TargetInfo(request)

    ## call GTC_AO to compute Strehl ratio and Encircled Energy
	guide_star = GuideStar(float(request.POST.get("gs_magnitude")),float(request.POST.get("gs_separation")))

	aocor = GTC_AO(sky_conditions,guide_star,telescope_params.aperture)

	# Creates an object of type Calculator_Spect
	a=Calculator_Spect(target_info,sky_conditions,spectrograph_setup,selected_scale,telescope=telescope_params, instrument=frida_setup, aocor = aocor)

	static_response, throughput = a.get_static_response()
	wave_array = a.grating_info.wave_array()

	atrans = a.get_atrans()
	grating_effic = a.get_grating_effic()
	sky_rad = a.get_sky_rad()

	context= {'debug_values':debug_values, "static_response":static_response, "throughput":throughput, "wave_array":wave_array, 'atrans': atrans, 'grating_effic':grating_effic, 'sky_rad':sky_rad}

	phot_zp = define_photzp()
	photons_obj = target_info.photons_wave(wave_array)
	photunit=u.photon/(u.s * u.cm**2 * u.angstrom)
	context['photons_obj'] = photons_obj.to(photunit)

	detector = frida_setup.detector
	phi_obj_total = photons_obj * atrans * telescope_params.area * telescope_params.reflectivity \
						 * grating_effic * throughput

	electron_per_s_unit=u.electron/(u.s * u.angstrom)
	exptime = 300 * u.s	## FIXME. What is this param in web?

	## compute sky emission and atmospheric absorption  (done within compute_sky_photonrate_spect)
	area_reference = np.pi*(50 * u.marcsec)**2
	phi_sky_apert = sky_rad * telescope_params.area * telescope_params.reflectivity * grating_effic * throughput * area_reference

	## compute S/N
	dit = float(request.POST.get('DIT_exp',1)) * u.s
	nexp = int(request.POST.get('N_exp',1)) 
	texp = nexp * dit
	aperture = {"EE":0.5,"Npix":27}		##FIXME. What is this param in web?
	phi_obj_apert = phi_obj_total * aperture["EE"]
	dwave = a.grating_info.delt_wave
	shot_noise2 = (texp * dwave * (phi_obj_apert + phi_sky_apert).to(electron_per_s_unit)).value
	dark_noise2 = (texp * detector['darkc'] * aperture['Npix']).value
	read_noise2 = (nexp * aperture['Npix'] * detector['ron'] ** 2).value
	noise = np.sqrt(shot_noise2 + dark_noise2 + read_noise2)
	signal = texp * dwave * phi_obj_apert
	snr = signal.to("electron").value / noise
	context["snr"] = snr
	debug_values.append("SNR => %s" % str(snr))


	debug_values = debug_values + target_info.debug()
	context['debug_values'] = debug_values

	return render (request, "calculate_ifs.html", context)


#	lambda_eff = a.wave_array[1000:1010]
#	strehl= aocor.compute_strehl(lambda_eff,sky_conditions['airmass'])
#	psf = aocor.compute_psf(lambda_eff,strehl['StrehlR'])
#	## the aperture in the IFS case is box-like, one side corresponds to the "slit-width" and the other
#	## side is along the radial profile
#	fcore = 1.5 # radius of aperture as a factor of FWHM of the PSF core
#	pixscale = frida_setup.pixscale # in arcseconds
#	ee_aperture=aocor.compute_ee(psf,pixscale,fcore=fcore,spaxel=1)
#
#	dit = float(request.POST.get('DIT_exp'))
#	Nexp = int(request.POST.get('N_exp'))
#	Nexp_min = 1
#	Nexp_max = 100.
#	Nexp_step = (Nexp_max-Nexp_min)/30
#	(texp_seq,snr_seq) = a.signal_noise_texp_spect(Nexp_min,Nexp_max,Nexp_step,dit,ee_aperture)
#	texp_seq = np.linspace(0,texp_seq, len(snr_seq))
#
#	required_sn = snr_seq[-1]	#ndit = a.texp_signal_noise_spect(required_sn,dit,aperture)
#	#ndit = a.texp_signal_noise_spect(required_sn,dit,aperture)
#	print (" Required S/N ",required_sn)
#	#print ("Time required to reach S/N ",ndit,ndit*dit)
#
#	################################
#	# Graphical output
#	# exposure time vs SNR
#	len_texp_seq = len(texp_seq)
#	m4 = np.zeros((len_texp_seq,2))
#	print ("texp, snr ",texp_seq[len_texp_seq-1],snr_seq[len_texp_seq-1])
#	m4[:,0] = texp_seq
#	m4[:,1] = snr_seq
#	m4 = str(m4.tolist())
#
#	# wavelength vs SED above atmosphere
##	nele = len(obs_filter.wave)
##	m2 = np.zeros((nele,2))
##	m2[:,0] = obs_filter.wave
##	m2[:,1] = target_info.flambda_wave(obs_filter.wave)/1.e-16
##	m2 = str(m2.tolist())
##
#
#	mag_target = float(request.POST.get('spatial_integrated_brightness'))
#	debug_values["FWHM_core"] = psf['FWHM_core']
#
#	context = {
#					'graph_title':'IFS Mode',
#					'sed_flambda':m4,
#					'signal_noise':m4,
#					'lambda_center': central_wave_grating,
#					'Object_magnitude': mag_target,
#					'Input_Band': target_info.Band,
#					'strehl_ratio':strehl['StrehlR'],
#					'encircled_energy': ee_aperture['EE'],
#					'Aperture_radius': ee_aperture['Radius'],
#					'Pixscale': pixscale,
#					'AreaNpix': ee_aperture['Npix'],
#					'flux_obj': None,
#					#'detected_photons_from_source': str("%9.2E"% a.phi_obj_total.value),
#					#'detected_photons_from_sky_sqarcsec' : str("%9.2E"% a.phi_sky_sqarc.value),
#					#'atmospheric_transmission': a.atmostrans[10],
#					'efficiency': a.grating_effic,
#					'debug_values': debug_values,
#					}
#	return render(request, 'result.html', context)
#
#

