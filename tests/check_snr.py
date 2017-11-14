import numpy as np
import matplotlib.pyplot as plt
from frida.set_config import *
from frida.gtcao import *
import scipy.constants as pconst


# telescope 
teldiam1=10.4  # in meter
teldiam2=1.
telarea=3.14*(teldiam1*teldiam1-teldiam2*teldiam2)/4.
refl_primary=0.9
refl_second=0.9
refl_tertia=0.9
effic_tel=refl_primary*refl_second*refl_tertia 

# Atmospheric conditions - fried parameter
#r0=10   
fwhm_seeing=0.9
lambda_seeing=0.55
airmass=1.2

## guide star characteristics
mag_gs=13.
sep_gs=0.

## instrument characteristics
filter='J'
param_filter = get_filter(filter)
lambda_obs = param_filter['lamb_eff']

pixscale = 0.01 # in arcseconds

param_instru = par_instru(filter)
param_atm = par_atm(filter)
param_gtcao = par_gtcao(filter)

## param_instru:
##   qe  - Quatum efficiency
##   effic_relay - Efficiency of relay optics
##   effic_came  - Efficiency of camera optics
effic_instr=param_instru['qe']*param_instru['effic_relay']*param_instru['effic_cam']*param_filter['trans']
effic_gtcao = param_gtcao['throughput']

## detector setup
dit = 10.
nexp = 10

fvega=flux_vega(filter)


## compute strehl and ee
strehl = compute_strehl(filter,fwhm_seeing,mag_gs,sep_gs,airmass)
ee_point = compute_ee(strehl['StrehlR'],fwhm_seeing,filter)
radaper = 1.5*ee_point['FWHM_core']
area_ref = pconst.pi*radaper*radaper  # arcsec^2
npix =  area_ref / pixscale / pixscale


## source characteristics
mag_obj=14.5
flux_obj = fvega['flux0']*10.**(-0.4*mag_obj)


# Vega zero point
#          B     V      R      I     J      H      Ks 
#lamb     0.438 0.545  0.641  0.798 1.22   1.63   2.19  microns  
#f_lamb   632  363.1  217.7  112.6  31.47  11.38 3.96  x10^-11 erg cm-2 s-1 A-1 
# J @ 1.25mic - 8.51 e W/m2/mic (log)

phi_obj_tel=telarea*effic_tel*param_atm['trans']*flux_obj/param_filter['energy_phot']*param_filter['bandwidth']
phi_obj_frida=phi_obj_tel*effic_gtcao*effic_instr



## sky emission
mag_sky=param_atm['mag_sqarc']
flux_sky = fvega['flux0']*10.**(-0.4*mag_sky)
phi_sky_tel=telarea*effic_tel*param_atm['trans']*flux_sky/param_filter['energy_phot']*param_filter['bandwidth']
phi_sky_frida=phi_sky_tel*effic_gtcao*effic_instr*area_ref

## compute S/N ratio
noise2_obj = dit*phi_obj_frida
noise2_sky = dit*phi_sky_frida
noise2_read = npix*param_instru['ron']**2
noise2_dark = dit*param_instru['dark']
snr = nexp*dit*phi_obj_frida/np.sqrt(nexp*(noise2_obj+noise2_sky+noise2_read+noise2_dark))

print(snr)
print(pconst.h*pconst.c)
print(param_atm['trans'])
print(flux_sky)
print(flux_obj)
print(lambda_obs)
print(param_filter)


 
