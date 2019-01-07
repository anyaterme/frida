from astropy.io import fits
import astropy.units as u

def add_keywords_target(hdr,target_info):
    hdr.set('InputMag',target_info.Magnitude,'Input target magnitude')
    hdr.set('MagSyst',target_info.MagSystem,'Magnitude system [AB/Vega]')
    hdr.set('SEDType',target_info.energy_type,'Flux distribution type')
    hdr.set('MagBand',target_info.Band,'Input mag. band')
    return hdr

def add_keywords_sky(hdr,sky_conditions):
    '''
    '''
    hdr.set('Seeing',sky_conditions['seeing'].value,'Seeing at reference wavelength')
    hdr.set('Airmass',sky_conditions['airmass'],'Airmass ')
    return hdr

def add_keywords_gtcao(hdr,guide_star):
    '''
    '''
    hdr.set('mag_GS',guide_star.Magnitude,'Magnitude of reference star')
    hdr.set('sep_GS',guide_star.Separation.value,'Axis distance of reference star')
    return hdr

def add_keywords_instru_ima(hdr,instrument_setup,obs_filter):
    '''
    '''
    hdr.set('ObsFilter',obs_filter.info['Name'],'Instrument Filter')
    hdr.set('pixscale',instrument_setup.pixscale.value,'Pixel scale [marcs]')
    return hdr

def add_keywords_instru_ifs(hdr,instrument_setup,spectrograph_setup):
    '''
    '''
    hdr.set('Grating',spectrograph_setup['Grating'],'Selected Grating')
    hdr.set('CentWave',spectrograph_setup['Central_wave'].value,'Mid Wavelength ['+spectrograph_setup['Central_wave'].unit.name+']')
    hdr.set('pixscale',instrument_setup.pixscale.value,'Pixel scale [mas]')
    return hdr

def write_fits_image(filename,data,target_info,sky_conditions,guide_star,instrument_setup,obs_filter):
    '''
    '''
    hdu = fits.PrimaryHDU(data)
    hdr=hdu.header
    hdr=add_keywords_target(hdr,target_info)
    hdr=add_keywords_sky(hdr,sky_conditions)
    hdr=add_keywords_gtcao(hdr,guide_star)
    hdr=add_keywords_instru_ima(hdr,instrument_setup,obs_filter)
    hdu.writeto(filename)

def write_fits_cube(filename,data,target_info,sky_conditions,guide_star,instrument_setup,spectrograph_setup,InitWave,deltaWave):
    '''
    '''
    hdu = fits.PrimaryHDU(data)
    hdr=hdu.header
    hdr.set('CTYPE1','x')
    hdr.set('CTYPE2','y')
    hdr.set('CTYPE3','WAVELENGTH')
    hdr.set('CUNIT3',InitWave.unit.name)
    hdr.set('CDELT1',instrument_setup.pixscale.value)
    hdr.set('CDELT2',instrument_setup.pixscale.value)
    hdr.set('CDELT3',deltaWave.value)
    hdr.set('CRPIX3',1)
    hdr.set('CRVAL3',InitWave.value)
    hdr=add_keywords_target(hdr,target_info)
    hdr=add_keywords_sky(hdr,sky_conditions)
    hdr=add_keywords_gtcao(hdr,guide_star)
    hdr=add_keywords_instru_ifs(hdr,instrument_setup,spectrograph_setup)
    hdu.writeto(filename)
