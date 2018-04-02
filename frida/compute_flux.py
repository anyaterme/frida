import numpy as np
import astropy.constants as const
import astropy.units as u

def energy_photon(wave):
	"""
        Compute the energy of a photon at a given wavelenght [in microns]
        	 :param wave:
	      :return: energy of a photon, units are energy / photon
	"""

	return const.h * const.c / wave / u.photon

def bbody(lamb, T):
    """ Blackbody as a function of wavelength ([lamb]=micron) and temperature ([T]=K).
        [h] = J s
        [k] = J K^-1
        [c] = m s^-1
    returns units of J s^-1 m^-3 - J s^-1 m^-2 * 10^-6 micron^-1    /s/cm^2/cm/Steradian
    """
    if not hasattr(lamb,'unit'):
        lamb = lamb * u.micron
    if not hasattr(T,'unit'):
        T = T * u.K
    lamb_m = lamb.to('m') # convert microns to metres
    #norm_J = 1.67e-16 ## normalization to have mag_Vega=15 at J with T=9600. which is Teff of A0
    norm_J = phot_zp['J']['l10_flamb'].physical ## normalization to have mag_Vega=15 at J with T=9600. which is Teff of A0
    exponent=(const.h*const.c / (lamb.si*const.k_B*T)).si
    flamb = norm_J * 2*const.h*const.c**2 / (lamb_m**5 * (np.exp(exponent) - 1))
    return flamb

def powl(lamb,alpha,band_norm='J'):
    """ Power law as a function of wavelength ([lamb]=micron) and alpha index .
        F_nu ~ nu^-alpha  ==> F_lamb ~ lambda ^ (alpha-2)
    returns units of J s^-1 m^-3 - J s^-1 m^-2 * 10^-6 micron^-1    /s/cm^2/cm/Steradian
    """
    lamb_norm = phot_zp[band_norm]['efflamb']  # lambda_ref is set to 1.22 micron
    norm_J = phot_zp[band_norm]['l10_flamb'].physical ## normalization to have mag_Vega=15 at J with T=9600. which is Teff of A0
    flamb = norm_J * (lamb/lamb_norm)**(alpha-2)
    return flamb

def stellar_template(lamb,spectype):
    flamb = lamb
    return flamb

def define_photzp():
    """
      Photometric zero points, first value is the reference wavelength
         second value is filter width [microns]
         third value is log10(F_lambda[Vega]) [W m^-2 micron^-1] [W m^-2 micron^-1]=10x[erg cm^-2 s^-1 Angstrom^-1]
         Zero_points are taken from Colina et al.
         fourth value is the difference between magnitudes m_AB - m_Vega
         To operate with dex use log10_flambZP.physical (quantity in scale and units)
    """
    phot_zp ={}
    
    id = 'B'
    efflamb = 0.438 * u.micron
    bwidth = 0.090 * u.micron
    log10_flambZP = -7.1938 * u.dex(u.Watt / u.m**2 / u.micron) 
    mAB_to_mVega = -0.09
    zp = {'name':id,'efflamb':efflamb,'bwidth':bwidth,'l10_flamb':log10_flambZP,'mAB_to_Vega':mAB_to_mVega}
    phot_zp[id]=zp    

    id = 'V'
    efflamb = 0.545 * u.micron
    bwidth = 0.085 * u.micron
    log10_flambZP = -7.4353 * u.dex(u.Watt / u.m**2 / u.micron) 
    mAB_to_mVega = 0.02
    zp = {'name':id,'efflamb':efflamb,'bwidth':bwidth,'l10_flamb':log10_flambZP,'mAB_to_Vega':mAB_to_mVega}
    phot_zp[id]=zp    

    id = 'R'
    efflamb = 0.641 * u.micron
    bwidth = 0.138 * u.micron       # Change by Dani. We cann't have two values equals for interpolate
    log10_flambZP = -7.7167 * u.dex(u.Watt / u.m**2 / u.micron) 
    mAB_to_mVega = 0.21
    zp = {'name':id,'efflamb':efflamb,'bwidth':bwidth,'l10_flamb':log10_flambZP,'mAB_to_Vega':mAB_to_mVega}
    phot_zp[id]=zp    

    id = 'I'
    efflamb = 0.798 * u.micron
    bwidth = 0.15 * u.micron
    log10_flambZP = -8.0273 * u.dex(u.Watt / u.m**2 / u.micron) 
    mAB_to_mVega = 0.45
    zp = {'name':id,'efflamb':efflamb,'bwidth':bwidth,'l10_flamb':log10_flambZP,'mAB_to_Vega':mAB_to_mVega}
    phot_zp[id]=zp    

    id = 'J'
    efflamb = 1.22 * u.micron
    bwidth = 0.26 * u.micron
    log10_flambZP = -8.5157 * u.dex(u.Watt / u.m**2 / u.micron) 
    mAB_to_mVega = 0.91
    zp = {'name':id,'efflamb':efflamb,'bwidth':bwidth,'l10_flamb':log10_flambZP,'mAB_to_Vega':mAB_to_mVega}
    phot_zp[id]=zp    

    id = 'H'
    efflamb = 1.63 * u.micron
    bwidth = 0.29 * u.micron
    log10_flambZP = -8.9586 * u.dex(u.Watt / u.m**2 / u.micron) 
    mAB_to_mVega = 1.39
    zp = {'name':id,'efflamb':efflamb,'bwidth':bwidth,'l10_flamb':log10_flambZP,'mAB_to_Vega':mAB_to_mVega}
    phot_zp[id]=zp    

    id = 'K'
    efflamb = 2.19 * u.micron
    bwidth = 0.41 * u.micron
    log10_flambZP = -9.4179 * u.dex(u.Watt / u.m**2 / u.micron) 
    mAB_to_mVega = 1.85
    zp = {'name':id,'efflamb':efflamb,'bwidth':bwidth,'l10_flamb':log10_flambZP,'mAB_to_Vega':mAB_to_mVega}
    phot_zp[id]=zp    
    
    return phot_zp

    
##phot_zp = {
##    """  Photometric zero points, first value is the reference wavelength
##         second value is filter width [microns]
##         third value is log10(F_lambda[Vega]) [W m^-2 micron^-1] [W m^-2 micron^-1]=10x[erg cm^-2 s^-1 Angstrom^-1]
##         Zero_points are taken from Colina et al.
##         fourth value is the difference between magnitudes m_AB - m_Vega
##    """
##    "B":[0.438, 0.090, -7.1938,-0.09],
##    "V":[0.545, 0.085, -7.4353, 0.02],
##    "R":[0.641, 0.15, -7.7167, 0.21],
##    "I":[0.798, 0.15, -8.0273, 0.45],
##    "J":[1.22, 0.26, -8.5157, 0.91],
##    "H":[1.63, 0.29, -8.9586, 1.39],
##    "K":[2.19, 0.41, -9.4179, 1.85]
##
##}

def scale_to_vega(band_ref,mag_ref,normflux_ref):
  """
      INPUTS: the magnitude at a reference band
  """
  print(band_ref)
  print(phot_zp[band_ref]['l10_flamb'])
  print(0.4*mag_ref)
  flux_vega_ref = phot_zp[band_ref]['l10_flamb'].physical * 10.**(-0.4*mag_ref)
  scale_factor = flux_vega_ref / normflux_ref
  return scale_factor

phot_zp = define_photzp()
