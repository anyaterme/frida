"""
A view function, or view for short, is simply a Python function that takes a Web request and returns a Web response.
This response can be the HTML contents of a Web page, or a redirect, or a 404 error, or an XML document, or
an image . . . or anything, really.
This module reads information from the web interface [params1.html-information from the target, 
params2.html-atmospheric conditions and guide star, and params3.html-instrument setup) and redirects the output to calculate_ima.html and
calculate_ifs.html

"""
from django.shortcuts import render
from django.conf import settings
from django.core.files import File
from django.http import HttpResponse
import glob
import os
import datetime, random, string, time
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import traceback
import astropy.units as u

from mpl_toolkits.mplot3d import Axes3D

#from math import sqrt, pi, log10

#import frida.calculator
#import frida.gtcao
#import frida.set_config
from frida.set_config import *
from frida.compute_flux import define_photzp
from frida.calculator import SkyConditions, GuideStar, TargetInfo, Calculator_Image, Calculator_Spect
from frida.gtcao import *
#from django.http import JsonResponse

# Create your views here.
def debug(msg=None, variable=None):
    if msg is None and variable is None:
        print ("%s @ HERE" % sys._getframe(1).f_code.co_name)
    elif variable is None:
        print ("%s @ %s"  % (sys._getframe(1).f_code.co_name, msg))
    elif msg is None:
        print ("%s @ %s"  % (sys._getframe(1).f_code.co_name, variable))
    else:
        print ("%s @ %s : %s"  % (sys._getframe(1).f_code.co_name, msg, variable))

def index(request):
    import csv

    #fp.close()
    #list_filters = {}


    fp1 = open(os.path.join(settings.INCLUDES, 'filters.dat'))
    list_filters = csv.DictReader(filter(lambda row: row[0] != '#', fp1))
    list_filters_dict = []
    telescope = "GTC"
    telescope_params = Telescope(telescope)
    for each_filter in list_filters:
        if (os.path.exists(os.path.join(settings.FILTERS, each_filter["Transmission"]))):
            my_filter = {}
            my_filter["Code"] = each_filter["Code"]
            my_filter["Name"] = each_filter["Name"]
            my_filter["Lambda"] = each_filter["lambda_center"] * u.micron 
            my_filter["Difr_limit"] = ((my_filter["Lambda"] / telescope_params.aperture) * u.rad).to(u.marcsec)
            print (my_filter["Difr_limit"])
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
        my_grating["Difr_limit"] = ((my_grating["Central_Wave"] * u.AA / telescope_params.aperture) * u.rad).to(u.marcsec)
        list_gratings_dict.append(my_grating)
    #fp2.close()
    fp3 = open(os.path.join(settings.SED_LIBRARY, 'pickles','pickles_options.csv'))
    list_pickles= csv.DictReader(filter(lambda row: row[0] != '#', fp3))
    list_pickles_dict = []
    for each_pickles in list_pickles:
        my_pickles = {}
        my_pickles["Label"] = each_pickles["Label"]
        my_pickles["sed_file"] = each_pickles["sed_file"]
        list_pickles_dict.append(my_pickles)
    fp4 = open(os.path.join(settings.SED_LIBRARY, 'nonstellar','nonstellar_options.csv'))
    list_nonstellar= csv.DictReader(filter(lambda row: row[0] != '#', fp4))
    list_nonstellar_dict = []
    for each_nonstellar in list_nonstellar:
        my_nonstellar = {}
        my_nonstellar["Label"] = each_nonstellar["Label"]
        my_nonstellar["sed_file"] = each_nonstellar["sed_file"]
        list_nonstellar_dict.append(my_nonstellar)


    context = {'list_filters':list_filters_dict, 'list_gratings':list_gratings_dict, 'list_pickles':list_pickles_dict, 'list_nonstellar':list_nonstellar_dict}
    print (os.path.abspath(os.path.dirname(__file__)))
    files = glob.glob("%s/*.png" % os.path.abspath(settings.MEDIA_ROOT))
    files.sort(key=os.path.getmtime)
    for f in files:
        try:
            if (time.time() - os.path.getmtime(f) >= 24*3600):
                os.remove(f)
        except Exception as e:
            print ("ERROR!!! File [%s] don't remove. :> %s" % (f, e))
    files = glob.glob("%s/*.fits" % os.path.abspath(settings.MEDIA_ROOT))
    files.sort(key=os.path.getmtime)
    for f in files:
        try:
            if (time.time() - os.path.getmtime(f) >= 24*3600):
                os.remove(f)
        except Exception as e:
            print ("ERROR!!! File [%s] don't remove. :> %s" % (f, e))
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
    :return: dict sky_conditions{'seeing','wave_seeing','airmass','pwv'}
    """
    # cond = Conditions()
    sky_conditions = {'seeing': 0.9 * u.arcsec, 'wave_seeing': 0.5 * u.micron, \
                   'airmass': 1.2, 'pwv': 2.5 * u.mm}
    sky_conditions['seeing'] = float(request.POST.get("seeing", sky_conditions['seeing'])) * u.arcsec
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
    error_values = {}
    # Request parameters relative to astronomical target
    mag_target = float(request.POST.get('spatial_integrated_brightness'))
    band = request.POST.get("band_flux")
    #mag_system = request.POST.get('units_sib')
    mag_system = request.POST.get('mag_system')
    energy_type = request.POST.get('spectral_type')
    A_V = request.POST.get('extinction')
    
    error = False

    ## Select morphology of source
    source_morpho = request.POST.get('source_type')
    if (source_morpho == 'extended'):
        label_source_morpho = 'Extended source'
        extended = True         
    elif (source_morpho == 'point'):
        label_source_morpho = 'Point source'
        extended = False        

    # create a tuple with the information relative to the Spectral Energy Distribution
    debug_values["Spectral Distribution"] = "Checking...."
    sed = (energy_type,0)
    if (energy_type == 'black_body'):
        try:
            temperature= float(request.POST.get('temperature'))
            sed = (energy_type,temperature)
            label_energy_type = 'Black Body, T=%s K' % temperature
            debug_values['Spectral Distribution'] = label_energy_type
        except:
            error_values['Spectral Distribution'] = 'Black Body, but you has indicated wrong temperature'
            error = True

    elif (energy_type == 'power_law'):
        try:
            pl_index = float(request.POST.get('pl_index'))
            sed = (energy_type,pl_index)
            label_energy_type = 'Power Law, lambda^%s' % pl_index
            debug_values['Spectral Distribution'] = label_energy_type
        except Exception as e:
            error_values['Spectral Distribution'] = 'Power Law, but you has indicated wrong lambda (%s)' % e
            error = True
    elif (energy_type == 'stellar_template'):
        try:
            templatename = request.POST.get('star_type')
            sed = ('pickles', templatename)
            label_energy_type = 'Pickles, template: %s' % templatename
            debug_values['Spectral Distribution'] = label_energy_type
        except:
            error_values['Spectral Distribution'] = label_energy_type
            error = True
    elif (energy_type == 'non-stellar-lib'):
        try:
            templatename = request.POST.get('st_non_stellar')
            sed = ('nonstellar', templatename)
            label_energy_type = 'Non-Stellar, template: %s' % templatename
            debug_values['Spectral Distribution'] = label_energy_type
        except:
            error_values['Spectral Distribution'] = label_energy_type
            error = True
    elif (energy_type == 'single_emission'):
        try:
            line_central_wave = request.POST.get('st_wavelength')
            line_flux = request.POST.get('st_flux')
            line_flux_units = request.POST.get('st_flux_units')
            if (line_flux_units == 'cgs'):
                line_flux = line_flux * u.erg / u.cm**2 / u.s
            elif (line_flux_units == 'si'):
                line_flux = line_flux * u.watt / u.m**2 / u.s
                  
            line_velocity = request.POST.get('st_velocity') * u.km / u.s
            sed = ('emission_line',line_central_wave,line_flux,line_velocity)
            label_energy_type = 'Emission Line at ' % templatename
            debug_values['Spectral Distribution'] = label_energy_type
        except:
            error_values['Spectral Distribution'] = label_energy_type
            error = True

    waveshift = None
    try:
        if (request.POST.get("velocity_type",'redshift') == 'redshift'):
            try:
                waveshift = ('redshift', 0.)
                my_value= request.POST.get('redshift_value', 0)
                if my_value == "":
                    my_value = 0
                waveshift = ('redshift', float(my_value))
                debug_values['Velocity Type'] = "Redshift, Z=%d" % float(my_value)
            except Exception as e:
                error_values['Velocity Type'] = "Redshift, but you has indicated wrong Z (%s)" % e
                error = True
        else:
            try:
                waveshift = ('radial_velocity', 0. * u.km / u.s )
                my_value = request.POST.get('radial_value', 0)
                if my_value == "":
                    my_value = 0
                waveshift = ('radial_velocity', float(my_value) * u.km / u.s)
                debug_values['Velocity Type'] = "Radial Velocity, v=%d [Km/s]" % float(my_value)
            except:
                error_values['Velocity Type'] = "Redshift, but you has indicated wrong velocity"
                error = True
    except:
        error_values['Velocity Type'] = "You must indicate velocity type"
        error = True
        pass

    # creates an object of type TargetInfo, it provides method to compute scaled f-lambda at any wavelength
    target_info = TargetInfo(mag_target,band,mag_system,sed,waveshift,extended=extended)
    target_info.error = error
    target_info.messages = debug_values
    target_info.error_messages = error_values
    target_info.energy_type = label_energy_type
    target_info.source_type= label_source_morpho
    
    return target_info

def calculate_ima(request):
    """
    Compute an object with the information to  about target from html templates. It returns an object of class TargetInfo
    :param request: input from html
    :return: object TargetInfo
    """
    debug_values = {}
    telescope = "GTC"
    telescope_params = Telescope(telescope)
    sky_conditions = get_SkyConditions(request)
    guide_star = GuideStar(float(request.POST.get("gs_magnitude")),\
                        float(request.POST.get("gs_separation"))*u.arcsec)
    selected_scale = request.POST.get('scale','fine_scale')
    #frida_setup = Instrument_static(request.POST.get('scale','fine_scale'),path=settings.INCLUDES)
    # Filter
    selected_filter = request.POST.get('filter')
    msg_err = {}
    try:
        obs_filter = Filter(selected_filter,path_list=settings.INCLUDES,path_filters=settings.FILTERS)
    except Exception as e:
        msg_err["ERROR GENERIC"] = ("%s" % e)
        return render(request, "error.html", {'msg_err':msg_err})
    lambda_eff=obs_filter.lambda_center()

    try:
        target_info = get_TargetInfo(request)
    except Exception as e:
        msg_err['ERROR GENERIC'] = "Please. review your TargetInfo parameters. %s" % e
        return render(request, "error.html", {'msg_err':msg_err})
    if (target_info.error):
        return render(request, "error.html", {'msg_err':target_info.error_messages})
    a=Calculator_Image(target_info,selected_filter,scale=selected_scale,\
                    telescope_name=telescope)
    static_response, throughput = a.get_static_response()

    print('@calculate_ima -> lambda_eff=',lambda_eff)    

    aocor = GTC_AO(sky_conditions,guide_star,telescope_params.aperture)
    strehl= aocor.compute_strehl(lambda_eff,sky_conditions['airmass'])
    psf = aocor.compute_psf(settings.PSF_MODEL,lambda_eff,strehl['StrehlR'])
    #fcore = 1.5 # radius of aperture as a factor of FWHM of the PSF core        FIX ME Where is the parameter?
    a.debug_values['Strehl='] = strehl['StrehlR'] 
    a.debug_values['FWHM_core='] = psf['FWHM_core'] 
    a.debug_values['FWHM_seeing='] = psf['FWHM_halo'] 
    pixscale = a.pixscale # in arcseconds

    ## Calculating output for Point source Extended source'
    fcore = float(request.POST.get('fcore', '1.5'))
    a.debug_values['fcore='] = fcore
    ## the aperture will be modified according to the geometry of the target
    ## if extended it will as reference area 1 pixel, fcore will be the
    aperture=aocor.compute_ee(psf,pixscale,fcore=fcore,\
                               source_type=target_info.source_type)

    if (target_info.source_type == "Point source"):
        fcore_seq = np.logspace(-0.1,2.,num=10)
        #fcore_seq = np.linspace(0.2,10,num=30)
        aperture_seq = aocor.compute_ee(psf,pixscale,fcore=fcore_seq,\
          source_type=target_info.source_type)
    
        for i in range(len(fcore_seq)):
            print("aperture_seq['Radius','EE']",aperture_seq['Radius'][i],\
                 aperture_seq['EE'][i]*100) 
   
    a.debug_values['Radius='] = aperture['Radius'] 
    a.debug_values['EE='] = aperture['EE'] 
    a.debug_values['EE-core='] = aperture['EE-core'] 
    a.debug_values['EE-halo='] = aperture['EE-halo'] 
    ## Calculating output for 'Extended source'

    dit = float(request.POST.get('DIT_exp','1')) * u.second
    Nexp = int(request.POST.get('N_exp','1'))
    dit = float(request.POST.get('DIT_exp')) * u.second
    signal_noise_dit = a.signal_noise_dit(dit,aperture)

    calc_method = request.POST.get('type_results')
    if (calc_method=="SN_ratio"):
        signal_noise= float(request.POST.get('signal_noise'))
        required_sn = signal_noise
        ndit = a.texp_signal_noise_img(required_sn,dit,aperture)
        Nexp = ndit
    else:
        Nexp = int(request.POST.get('N_exp'))

    Nexp_min = max(1,Nexp-2)
    Nexp_max = Nexp+2
    Nexp_step = (Nexp_max-Nexp_min)/5
    signal_noise_seq = a.signal_noise_texp_img(Nexp_min,Nexp_max,Nexp_step,dit,aperture)
    texp_seq = signal_noise_seq['texp']
    snr_seq = signal_noise_seq['SNR']
    aux = np.abs(texp_seq - Nexp*dit)
    signal_required_sn = signal_noise_seq['signal'][aux.argmin()]
    required_sn = snr_seq[aux.argmin()]

    '''
    Nexp = int(request.POST.get('N_exp'))
    Nexp_min = max(1,Nexp-10)
    Nexp_max = Nexp+20
    Nexp_step = (Nexp_max-Nexp_min)/30
    signal_noise_seq = a.signal_noise_texp_img(Nexp_min,Nexp_max,Nexp_step,dit,aperture)


    ## JAP - FIXME no esta claro que sea necesario 
    signal_noise= float(request.POST.get('signal_noise'))
    if (calc_method=="SN_ratio"):
        required_sn = signal_noise
    else:
        required_sn = snr_seq[np.where(texp_seq == Nexp * dit)]
    try:
        ndit = a.texp_signal_noise_img(required_sn,dit,aperture)
        if (calc_method == "SN_ratio"):
            Nexp = ndit[0]
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print (traceback.print_exc())
        print (exc_type, exc_obj, exc_tb.tb_frame.f_code.co_filename, exc_tb.tb_lineno)
        msg_err["ERROR GENERIC"]=("%s" % e)
        return render(request, "error.html", {'msg_err':msg_err})
    '''

    ## compute image using PSF information 
    limit_window = (settings.IMSIZE_IMG,settings.IMSIZE_IMG)
    print ("######", limit_window)
    if (settings.PSF_MODEL == "Airy+Gaussian"): 
       im_psf2d,im_xaxis,im_yaxis =buildim_psf_AiryGauss(psf,a.pixscale,Nx=limit_window[0],Ny=limit_window[1])
    elif (settings.PSF_MODEL == "2-Gaussians"):  
       im_psf2d,im_xaxis,im_yaxis =buildim_psf_2gauss(psf,a.pixscale,Nx=limit_window[0],Ny=limit_window[1])

    if hasattr(im_psf2d,'unit'):
         print('units im_psf2d ',im_psf2d.unit)
    else:
         print('no units im_psf2d ')

    ## now scale with the signal, psf is normalize to have area unity 
    if (target_info.source_type == "Point source"):
        print ("#########################################################################")
        im_signal_obj = a.phi_obj_total*im_psf2d * aperture['Area_pixel']
        max_signal_obj = np.max(im_signal_obj)
        print (im_signal_obj.unit)
        print (max_signal_obj.unit)
        print ( (a.phi_sky_sqarc*aperture['Area_pixel']).unit)
        max_signal_obj_sky = max_signal_obj + a.phi_sky_sqarc*aperture['Area_pixel']
    elif (target_info.source_type == "Extended source"):        
        max_signal_obj = a.phi_obj_total*aperture['Area_pixel']
        max_signal_obj_sky = max_signal_obj + a.phi_sky_sqarc*aperture['Area_pixel']

    max_signal_obj_sky_dit = dit * max_signal_obj_sky
    saturation_time = a.instrument.detector["welldepth"] / max_signal_obj_sky       

    print("max_signal_obj,obj_sky=",max_signal_obj.to(u.electron/u.second),\
       max_signal_obj_sky.to(u.electron/u.second))  
    
    ## create png
    context = {}
    name = None
    if (request.POST.get("2d_psf", "off") == "on") and (target_info.source_type == "Point source"):
        im_value = im_signal_obj.value
        lim_value = np.log10(im_value)
        print("Median Min max image",np.percentile(im_value,50),np.min(im_value),np.max(im_value))
        im_value /= np.max(im_value)        
        vmin = np.percentile(lim_value,15)
        vmax = np.percentile(lim_value,85)
        print("VMin Vmax image",vmin,vmax)
        print("Median Min max image",np.percentile(im_value,50),np.min(im_value),np.max(im_value))
        fig = plt.figure()
        ax = Axes3D(fig,elev=settings.IM3D_ELEV,azim=settings.IM3D_AZIM)
        # make the panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        xx, yy = np.meshgrid(im_xaxis,im_yaxis)        
        ax.plot_surface(xx,yy,im_value,rstride=1,cstride=1,cmap='hot')
        ax.set_zlim(-2,1)
        ax.contourf(xx,yy,im_value,zdir='z',offset=-2,cmap='hot',vmax=0.8,vmin=0) #, vmin=vmin,vmax=vmax)        
        #plt.imshow(im_value, vmin=vmin,vmax=vmax, cmap='hot')
        ##plt.imshow(lim_value), vmin=vmin,vmax=vmax, cmap='hot')
         ###ax.set_xlabel('Arcsec')
        ###ax.set_ylabel('Arcsec')
        name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
        print('file name=',name)
        f=open(os.path.join(settings.MEDIA_ROOT,'%s.png' % name), 'wb')
        plt.savefig(f, dit=2000)
        f.close()
        np.save(os.path.join(settings.MEDIA_ROOT,'%s.npy' % name), im_signal_obj.value)
        createfits(name)
        context['img_name'] = '%s' % name
        ## now scale with the signal, psf is normalize to have area unity 
    print('phi_sky_sqarc.unit',a.phi_sky_sqarc)
    print('Unit aperture[Area_pixel]=',aperture['Area_pixel'].unit)

    dit_pattern = float(request.POST.get('dit_pattern', '1.0'))
    
    context['Object_magnitude'] = target_info.Magnitude,
    context['dit_pattern'] = dit_pattern
    context['target_info'] = target_info
    context['filter_trans_wave'] = a.img_wave.to(u.AA)
    context['filter_trans'] = a.img_filter_trans 
    context['sky_conditions'] = sky_conditions
    context['guide_star'] = guide_star
    context['frida_setup'] = a.instrument
    context['Input_Band'] = target_info.Band
    context['strehl_ratio'] = strehl['StrehlR']
    context['encircled_energy'] = aperture["EE"]
    context['Aperture_radius'] = aperture['Radius']
    context['encircled_energy_seq'] = aperture_seq["EE"]
    context['Aperture_radius_seq'] = aperture_seq['Radius']
    context['pixscale'] = pixscale
    context['AreaNpix'] = aperture['Npix']
    context['static_response'] = static_response
    context['throughput'] = throughput
    context['throughput_lambda'] = throughput[(np.abs(a.img_wave-obs_filter.wave_median)).argmin()]
    context['atrans'] = a.atmostrans
    context['sky_rad'] = a.skyemission_photons
    context['global_effic'] = a.throughput 
    context['wave_array'] = a.img_wave.to(u.AA)
    context['texp_array'] = texp_seq
    context['snr'] = snr_seq
    context["obs_filter"] = obs_filter
    context['signal_req'] = signal_required_sn
    context['signal_req_per_dit'] = signal_required_sn / Nexp
    context['max_signal_obj_sky_dit'] = max_signal_obj_sky_dit 
    context['signal_noise_req'] = required_sn
    context['total_exposure_time'] = dit * Nexp
    context['dit'] = dit
    context['ndit'] = Nexp
    context['saturation_time'] = saturation_time.to(u.s)
    context['BLIP_time'] = a.blip_time
    context['background'] = a.phi_sky_sqarc*aperture['Area_pixel']*dit
    context['fwhm_core'] = psf['FWHM_core'].to(u.mas)
    context['darkc'] = a.instrument.detector['darkc'] * dit
    
    context['snr_graph'] = (request.POST.get("sn_as_exp_time", "off") == "on")
    context['EE_graph'] = (request.POST.get("enc_energy_aperture_radius", "off") == "on")
    context['phot_obj_trans'] = a.objphot_wave*a.atmostrans
    context['phot_obj'] = a.objphot_wave
    context['flambda_obj'] = a.flambda_wave
    print('signal_noise_dit:',signal_noise_dit['phot_obj'])
    print('a.phi_obj',a.objphot_wave)


    a.debug_values['throughput_lambda_index'] =(np.abs(obs_filter.wave-obs_filter.wave_median)).argmin()
    a.debug_values['Filter Wave'] = a.img_wave.to(u.AA)
    a.debug_values['Filter Trans Wave'] = a.img_filter_trans
    context['debug_values'] = a.debug_values

    return render(request, 'calculate_ima.html', context)

def createfits(filename):
    from astropy.io import fits
    from django.utils.encoding import smart_str
    data = np.load(os.path.join(settings.MEDIA_ROOT,'%s.npy' % (filename)))
    hdul = fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    hdul.append(fits.ImageHDU(data=data))
    hdul.writeto(os.path.join(settings.MEDIA_ROOT,'%s.fits' % (filename)))
    return None

def downloadfits(request, pngfile):
    from astropy.io import fits
    from django.utils.encoding import smart_str
    data = np.load(os.path.join(settings.MEDIA_ROOT,'%s' % (pngfile.replace('png','npy'))))
    hdul = fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    hdul.append(fits.ImageHDU(data=data))
    hdul.writeto(os.path.join(settings.MEDIA_ROOT,'%s' % (pngfile.replace('png','fits'))))
    response= HttpResponse(mimetype='application/force-download')
    response['Content-Disposition'] = 'attachment; filename=%s' % pngfile.replace('png','fits')
    response['X-Sendfile'] = os.path.join(settings.MEDIA_ROOT,'%s' % (pngfile.replace('png','fits')))
    return response



def calculate_ifs(request,telescope=settings.TELESCOPE):
    debug_values = ["test ==> test_content"]

    dit_pattern = float(request.POST.get('dit_pattern', "1."))
    telescope_params = Telescope(telescope)
    sky_conditions = get_SkyConditions(request)
    guide_star = GuideStar(float(request.POST.get("gs_magnitude")),\
                        float(request.POST.get("gs_separation"))*u.arcsec)

    selected_scale = request.POST.get('scale','fine_scale')
    
    #frida_setup = Instrument_static(scale=selected_scale,instid=settings.INSTRUMENT,path=settings.INCLUDES)

    selected_grating = request.POST.get('grating')
    central_wave_grating = float(request.POST.get('grating_value',None))*u.AA
    spectrograph_setup = {'Grating':selected_grating,'Central_wave':central_wave_grating}

    # creates an object of type TargetInfo, it provides method to compute scaled f-lambda at any wavelength
    target_info = get_TargetInfo(request)

    print('@calculate_ifs -> central_wave_grating=',central_wave_grating)    
    # Creates an object of type Calculator_Spect
    a=Calculator_Spect(target_info,spectrograph_setup,scale=selected_scale,\
            telescope_name=telescope)

    wave_array = a.wave_array
    atrans = a.atm_trans
    grating_effic = a.grating_effic
    sky_rad = a.sky_rad
    
    pixscale = a.pixscale

    ## call GTC_AO to compute Strehl ratio and Encircled Energy
    aocor = GTC_AO(sky_conditions,guide_star,telescope_params.aperture)
    strehl = aocor.compute_strehl(central_wave_grating,sky_conditions['airmass'])
    psf = aocor.compute_psf(settings.PSF_MODEL,central_wave_grating,strehl['StrehlR'])

    ## Calculating output for Point source Extended source'
    fcore = float(request.POST.get('fcore', '1.5'))
    a.debug_values['fcore='] = fcore
    ## the aperture will be modified according to the geometry of the target
    ## if extended it will as reference area 1 pixel, fcore will be the
    aperture=aocor.compute_ee(psf,pixscale,fcore=fcore,\
                source_type=target_info.source_type,spaxel=True)

    dit = float(request.POST.get('DIT_exp','10')) * u.second
    wave_ref = float(request.POST.get('ifs_lambda_ref')) * u.angstrom
    dwave_ref = 10. * u.angstrom  ## FIXME (should be obtained from WEB-interface)
    calc_method = request.POST.get('type_results')
    signal_noise= float(request.POST.get('signal_noise'))
    if (calc_method=="SN_ratio"):
        required_sn = signal_noise
        ndit = a.texp_signal_noise_spect(required_sn,dit,aperture,wave_ref,dwave_ref)
        print('ndit',ndit)
        Nexp = ndit
    else:
        Nexp = int(request.POST.get('N_exp','1'))


    '''
    calc_method = request.POST.get('type_results')
    signal_noise= float(request.POST.get('signal_noise'))
    required_sn = signal_noise
    dit = float(request.POST.get('DIT_exp','10')) * u.second
    Nexp = int(request.POST.get('N_exp','1'))
    signal_noise= float(request.POST.get('signal_noise'))
    '''
    ## compute S/N at reference wavelength at different exposure times

    Nexp_min = max(1,Nexp-2)
    Nexp_max = Nexp+2
    Nexp_step = (Nexp_max-Nexp_min)/5
    #signal_noise_seq = a.signal_noise_texp_seq_spect(Nexp_min,Nexp_max,Nexp_step,\
    #                        dit,aperture,wave_ref,dwave_ref)
    #texp_seq = signal_noise_seq['texp']
    #snr_seq = signal_noise_seq['SNR']
    #signal_seq = signal_noise_seq['signal']
    signal_noise_nexp = a.signal_noise_texp_spect(Nexp,dit,aperture)
    texp_seq = signal_noise_nexp['texp']
    snr_seq = signal_noise_nexp['SNR']
    signal_seq = signal_noise_nexp['signal']
    
    print("@calculate_ifs -> snr_seq",snr_seq)

    ## compute AO correction for the wave_array
    strehl_wave_array= aocor.compute_strehl(wave_array,sky_conditions['airmass'])
    psf_wave_array = aocor.compute_psf(settings.PSF_MODEL,wave_array,strehl['StrehlR'])

    ## compute image using PSF information 
    limit_window = (settings.IMSIZE_IFS,settings.IMSIZE_IFS)
    print ("######", limit_window)
    if (settings.PSF_MODEL == "Airy+Gaussian"): 
       im_psf_wave, cube_xaxis, cube_yaxis =buildcube_psf_AiryGauss(psf_wave_array,a.pixscale,
                    Nx=limit_window[0],Ny=limit_window[1])
    elif (settings.PSF_MODEL == "2-Gaussians"):  
       im_psf_wave, cube_xaxis, cube_yaxis =buildcube_psf_2gauss(psf_wave_array,a.pixscale,\
                    Nx=limit_window[0],Ny=limit_window[1])

    ## now scale with the signal, psf is normalize to have area unity 
    cube_signal_obj = np.zeros_like(im_psf_wave) * a.phi_obj_spect[0].unit
    print ("*********************************************************************************************")
    print (im_psf_wave.unit)
    print (a.phi_obj_spect[0].unit)
    if (target_info.source_type == "Point source"):
        for iwave in range(len(wave_array)):
            cube_signal_obj[iwave,:,:] = a.phi_obj_spect[iwave]*im_psf_wave[iwave,:,:]
        max_signal_obj = np.amax(cube_signal_obj)  
        print (max_signal_obj.unit)
        print((a.phi_sky_spect_sqarc*aperture['Area_pixel']).unit)
        max_signal_obj_sky = max_signal_obj + a.phi_sky_spect_sqarc*aperture['Area_pixel']
    elif (target_info.source_type == "Extended source"):        
        max_signal_obj = a.phi_obj_spect*aperture['Area_pixel']
        max_signal_obj_sky = max_signal_obj + a.phi_sky_spect_sqarc*aperture['Area_pixel']

    ## create png
    context = {}
    name = None
    if (request.POST.get("2d_psf", "off") == "on") and (target_info.source_type == "Point source"):
        indwave = 1000
        im_value = cube_signal_obj[indwave,:,:].value   
        vmin = np.percentile(im_value,25)
        vmax = np.percentile(im_value,85)
        fig = plt.figure()
        ax = Axes3D(fig,elev=settings.IM3D_ELEV,azim=settings.IM3D_AZIM)
        xx, yy = np.meshgrid(im_xaxis,im_yaxis)        
        ax.plot_surface(xx,yy,im_value,rstride=1,cstride=1,cmap='hot')
        ax.set_zlim(-3,2)
        ax.contourf(xx,yy,im_value,zdir='z',offset=-3,cmap='hot',vmin=vmin,vmax=vmax)        
        #plt.imshow(im_signal_obj.value, vmin=vmin,vmax=vmax, cmap='hot')
        ax.set_xlabel('Arcsec')
        ax.set_ylabel('Arcsec')
        name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
        print('file name=',name)
        f=open(os.path.join(settings.MEDIA_ROOT,'%s.png' % name), 'wb')
        #plt.savefig(f, dit=2000)
        fig.savefig(f,dit=2000)        
        f.close()
        context['img_name'] = '%s.png' % name
        ## now scale with the signal, psf is normalize to have area unity 


    context= {'debug_values':debug_values, "static_response":a.static_response, \
           "throughput":a.throughput, "wave_array":wave_array, 'atrans': atrans, \
           'grating_effic':grating_effic, 'sky_rad':sky_rad}

    photunit = u.electron/u.second/u.Angstrom

    lambda_ref = central_wave_grating
    context['guide_star'] = guide_star
    context['photons_obj'] = a.phi_obj_spect
    context['frida_setup'] = a.instrument
    context['target_info'] = target_info
    context['sky_conditions'] = sky_conditions
    context['signal_noise_req'] = signal_noise
    context['total_exposure_time'] = dit * Nexp
    context['Aperture_radius'] = aperture['Radius']
    context['encircled_energy'] = aperture["EE"]
    context['lambda_eff'] = lambda_ref
    context['throughput_lambda'] = a.throughput[(np.abs(wave_array-lambda_ref)).argmin()]
    context['grating'] = spectrograph_setup
    context['dit'] = dit
    context['ndit'] = Nexp
    context['strehl_ratio'] = strehl['StrehlR']
    context["snr"] = snr_seq
    context["signal"] = signal_seq
    debug('signal_seq', signal_seq)
    context['AreaNpix'] = aperture['Npix']
    context['dit_pattern'] = dit_pattern
    context['darkc'] = a.instrument.detector['darkc'] * dit

    debug_values = debug_values + target_info.debug()
    context['debug_values'] = debug_values

    return render (request, "calculate_ifs.html", context)


def docs_p1(request):
    context = {}
    return render (request, "doc-params1.html", context)

def docs_p2(request):
    context = {}
    return render (request, "doc-params2.html", context)

def docs_p3(request):
    context = {}
    return render (request, "doc-params3.html", context)
