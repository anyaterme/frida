# load needed modules
from frida.set_config import *
from frida.calculator import TargetInfo
import astropy.units as u

import numpy as np

from frida.compute_flux import *
import matplotlib.pyplot as plt


# Telescope & instrument setup
data_path = settings.INCLUDES 

phot_zp = define_photzp()

wave_array=np.arange(8900,25000.,5.)*u.angstrom

## target info
print("Set target parameters")
mag_target = 17.
band = 'H'
mag_system='Vega'
redshift=1.5
waveshift=('redshift',redshift)
print('band:',band,' ZP=',phot_zp[band],' bwidth=',phot_zp[band]['bwidth'])

print ("  Object mag=",mag_target," band ",band," redshift",redshift)


print(" .................. ")

print (" SED BlackBody T=5000 K, z=1.5")

sed=('black_body',5000.)

# sed=('power_law',0.)
target_info = TargetInfo(mag_target,band,mag_system,sed,waveshift)

sed_wave = target_info.sed_wave
print(' sed_wave ',sed_wave)

sed_flambda = target_info.sed_flambda
print(' sed_flambda ',sed_flambda)

flambunit=u.erg/(u.s * u.cm**2 * u.angstrom) 

flamb_obj =target_info.flambda_wave(wave_array) 
print(' flamb_obj ',flamb_obj.to(flambunit))


plt.plot(sed_wave.to('angstrom'),sed_flambda.to(flambunit),'.r')
plt.plot(wave_array.to('angstrom'),flamb_obj.to(flambunit),'xg')
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('r'+flambunit.to_string('latex'),size=2)
plt.savefig('tests/BB5000Kz1p5_flambda.png')
plt.close()


photons_obj=target_info.photons_wave(wave_array)

photunit=u.photon/(u.s * u.cm**2 * u.angstrom)
print(' photons_obj ',photons_obj[0:10].to(photunit))

plt.plot(wave_array.to('angstrom'),photons_obj.to(photunit),'.r')
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('r'+photunit.to_string('latex'),size=20)
plt.savefig('tests/BB5000Kz1p5_photons.png')
plt.close()


print(" .................. ")

print (" SED power law alpha=1.5")

sed=('power_law',1.5)

# sed=('power_law',0.)
target_info = TargetInfo(mag_target,band,mag_system,sed,waveshift)

sed_wave = target_info.sed_wave
print(' sed_wave ',sed_wave)

sed_flambda = target_info.sed_flambda
print(' sed_flambda ',sed_flambda)

flamb_obj =target_info.flambda_wave(wave_array) 

flambunit=u.erg/(u.s * u.cm**2 * u.angstrom) 
print(' flamb_obj ',flamb_obj.to(flambunit))


plt.plot(sed_wave.to('angstrom'),sed_flambda.to(flambunit),'.r')
plt.plot(wave_array.to('angstrom'),flamb_obj.to(flambunit),'xg')
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('r'+flambunit.to_string('latex'),size=2)
plt.savefig('tests/PwL1p5z1p5_flambda.png')
plt.close()



photons_obj=target_info.photons_wave(wave_array)

photunit=u.photon/(u.s * u.cm**2 * u.angstrom)
print(' photons_obj ',photons_obj.to(photunit))

plt.plot(wave_array.to('angstrom'),photons_obj.to(photunit),'.r')
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('r'+photunit.to_string('latex'),size=20)
plt.savefig('tests/PwL1p5z1p5_photons.png')
plt.close()


'''
photunit=u.photon/(u.s * u.cm**2 * u.angstrom)
print ("photons_obj=",photons_obj[1000:1010].to(photunit))
plt.plot(wave_array.to('angstrom'),photons_obj.to(photunit))
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
#plt.ylabel(r'$\rm{\ L_{\lambda} [{photons\ s^{-1} }]}$',size=20)
plt.ylabel('r'+photunit.to_string('latex'),size=20)
plt.savefig('tests/BB5000Kz1p5_photons.png')
plt.close()
'''

print(" .................. ")

print (" SED Pickles K0V (Teff=5240 K) ")

sed=('pickles','ukk0v.dat')

# sed=('power_law',0.)
target_info = TargetInfo(mag_target,band,mag_system,sed,waveshift)

sed_wave = target_info.sed_wave
print('unit sed_wave ',sed_wave)

sed_flambda = target_info.sed_flambda
print(' sed_flambda ',sed_flambda)

flamb_obj =target_info.flambda_wave(wave_array) 

flambunit=u.erg/(u.s * u.cm**2 * u.angstrom) 
print(' flamb_obj ',flamb_obj[0:10].to(flambunit))


plt.plot(sed_wave.to('angstrom'),sed_flambda.to(flambunit),'.r')
plt.plot(wave_array.to('angstrom'),flamb_obj.to(flambunit),'xg')
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('r'+flambunit.to_string('latex'),size=2)
plt.savefig('tests/ukk0vz1p5_flambda.png')
plt.close()


photons_obj=target_info.photons_wave(wave_array)

photunit=u.photon/(u.s * u.cm**2 * u.angstrom)
print(' photons_obj ',photons_obj[0:10].to(photunit))


plt.plot(wave_array.to('angstrom'),photons_obj.to(photunit))
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('r'+photunit.to_string('latex'),size=20)
plt.savefig('tests/ukk0vz1p5_photons.png')
plt.close()


print(" .................. ")

print (" SED QSO2 ")

sed=('nonstellar','qso2.dat')

# sed=('power_law',0.)
target_info = TargetInfo(mag_target,band,mag_system,sed,waveshift)

sed_wave = target_info.sed_wave
print('unit sed_wave ',sed_wave)

sed_flambda = target_info.sed_flambda
print(' sed_flambda ',sed_flambda)

flamb_obj =target_info.flambda_wave(wave_array) 

flambunit=u.erg/(u.s * u.cm**2 * u.angstrom) 
print(' flamb_obj ',flamb_obj[0:10].to(flambunit))


plt.plot(sed_wave.to('angstrom'),sed_flambda.to(flambunit),'.r')
plt.plot(wave_array.to('angstrom'),flamb_obj.to(flambunit),'xg')
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('r'+flambunit.to_string('latex'),size=2)
plt.savefig('tests/qso2z1p5_flambda.png')
plt.close()

photons_obj=target_info.photons_wave(wave_array)

photunit=u.photon/(u.s * u.cm**2 * u.angstrom)
print(' photons_obj ',photons_obj.to(photunit))


plt.plot(wave_array.to('angstrom'),photons_obj.to(photunit))
plt.xlabel(r'$\lambda \rm{[\AA ]}$')
plt.ylabel('r'+photunit.to_string('latex'),size=20)
plt.savefig('tests/qso2z1p5_photons.png')
plt.close()

