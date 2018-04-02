import numpy as np
from math import ceil
import os,re
import csv
import astropy.units as u
import scipy.integrate


def read_sed_pickles_files(fname_sed,path_sed='sed_library/pickles'):
	pickles=np.loadtxt(os.path.join(path_sed, fname_sed))
	#print(lines_sed.fieldnames[1])
	unitx=u.Unit('angstrom')
	unity=u.Unit('erg s-1 cm-2 angstrom-1')
	wave = np.array(pickles[:,0])
	## original SED are normalized to flux at 5500AA, here we arbitrarily scale to about 
	##  V~15, flambda=3.5e-15 erg/cm2/s/A
	flambda= np.array(pickles[:,1]) * 3.5e-15
	sed = {'wave':wave*unitx,'flambda':flambda*unity}
	return sed 

def read_sed_nonstellar_files(fname_sed,path_sed='sed_library/nonstellar'):
	fp = open(os.path.join(path_sed, fname_sed))
	lines_sed = csv.DictReader(filter(lambda row: row[0] != '#', fp),delimiter=" ",\
	  skipinitialspace=True)
	print(lines_sed.fieldnames[1])
	unitx=u.Unit(re.search('\[(.+?)\]',lines_sed.fieldnames[0]).group(1))
	unity=u.Unit(re.search('\[(.+?)\]',lines_sed.fieldnames[1]).group(1))
	print("unitx,y",unitx,unity)
	sed_wave = []
	sed_value = []
	for line in lines_sed:
		sed_wave.append(float(line[lines_sed.fieldnames[0]]))
		sed_value.append(float(line[lines_sed.fieldnames[1]]))
	## original SED are normalized to flux at 5500AA, here we arbitrarily scale to about 
	##  V~15, flambda=3.5e-15 erg/cm2/s/A
	flambda= np.array(sed_value) * 3.5e-15 	
	sed = {'wave':np.array(sed_wave)*unitx,'flambda':flambda*unity}
	fp.close()
	return sed


def read_grating_files(fname_gratings,path_gratings='gratings'):
	fp = open(os.path.join(path_gratings, fname_gratings))
	lines_fp = csv.DictReader(filter(lambda row: row[0] != '#', fp),delimiter=",")
	unitx=u.Unit(re.search('\[(.+?)\]',lines_fp.fieldnames[0]).group(1))
	unity=u.Unit(re.search('\[(.+?)\]',lines_fp.fieldnames[1]).group(1))
	gr_wave = []
	gr_value = []
	for line in lines_fp:
		gr_wave.append(float(line[lines_fp.fieldnames[0]]))
		gr_value.append(float(line[lines_fp.fieldnames[1]]))
	if (unity == '%'):	
		gr_value = np.array(gr_value) / 100.
	
	grating = {'wave':np.array(gr_wave)*unitx,'effic':np.array(gr_value)}
	fp.close()
	return grating


def convolve2gauss(wave,response,gauss_width):

	# get the number of pixels
	delt_wave = wave[1]-wave[0]
	gauss_width_units = ceil(gauss_width / delt_wave)
	nwidth  = 4
	xpulso = np.arange(0,nwidth*gauss_width_units)
	xcenter = nwidth*gauss_width_units/2.
	pulsogauss = np.exp(-(xpulso-xcenter)**2/gauss_width_units**2)
	pulsogauss = pulsogauss / pulsogauss.sum()
	conv_response = np.convolve(response,pulsogauss,'same')
	return conv_response


def convolve2pulse(wave,response,pulse_width):

	# get the number of pixels
	delt_wave = wave[1]-wave[0]
	pulse_width_units = int(ceil(pulse_width / delt_wave))
	pulsodelta = np.zeros(3*pulse_width_units)
	for i in range(pulse_width_units):
		pulsodelta[i+pulse_width_units] = 1
	pulsodelta = pulsodelta / pulsodelta.sum()
	conv_response = np.convolve(response,pulsodelta,'same')
	return conv_response


def interpolate(wvl, specfl, new_wvl, unity='None'):
	#from scipy.interpolate import InterpolatedUnivariateSpline as inter
#	from numpy import interp as inter
#	from scipy import polyval,polyfit
#	try:
#		
#		specfl = specfl[0::len(specfl)/len(wvl)]
#		specfl = specfl[0:len(wvl)]				
#		if (type(wvl) != type (new_wvl)):
#			if (hasattr(wvl, "unit")):
#				new_wvl = new_wvl * wvl.unit
#			else:
#				wvl = wvl * new_wvl.unit
#		print (np.min(wvl), np.min(new_wvl), np.max(new_wvl), np.max(wvl))
#		if (np.min(wvl) <= np.min(new_wvl)) and (np.max(new_wvl) <= np.max(wvl)):
#			if (unity == 'perone'):
#					result=inter(new_wvl, wvl, specfl).clip(0,1)
#			elif (unity == 'percent'):
#					result=inter(new_wvl, wvl, specfl).clip(0,100)
#			else:
#					result=inter(new_wvl, wvl, specfl)
#			return result
#		else:
#			print 3
#			print (len(wvl), len(specfl))
#			inter_func=np.polyfit(wvl, specfl, 1)
#			f = np.poly1d(inter_func)
#			return f(new_wvl)
	try:
		z = np.polyfit(wvl, specfl,1 )
		f = np.poly1d(z)
		return (f(new_wvl))
	except Exception as e:
		print ("ERROR in interpolate: ", e)
		return None

def interpolate_old(wvl,specfl,new_wvl,unity='None'):
	#
	#	   Interpolation to a given wvl array. Returned values
	#	   are clipped to 0, and if it is a transmission curve,
	#	   also to ones.
	#
	#	   Order=1 is used, as larger orders produce
	#	   weird effects around corners.
	#
	from scipy.interpolate import InterpolatedUnivariateSpline as inter
	from scipy import polyval,polyfit
	print (len(wvl), len(specfl))
	print (np.max(wvl), np.min(wvl))
	print (np.max(new_wvl), np.min(new_wvl))
	#
	inter_func=inter(wvl,specfl,k=1)
	if (unity == 'perone'):
			result=inter_func(new_wvl).clip(0,1)
	elif (unity == 'percent'):
			result=inter_func(new_wvl).clip(0,100)
	else:
			result=inter_func(new_wvl)
	#
	#	   To avoid weird extrapolations, a lineal fit to the
	#	   extreme points of the original curve is used
	#
	working_unit = None
	if hasattr(wvl, 'unit'):
		working_unit = wvl.unit
		wvl = np.array([item.value for item in wvl])
	ind_i=np.where(new_wvl <= wvl.min())[0]
	ind_f=np.where(new_wvl >= wvl.max())[0]
	nel=int(np.clip(0.1*len(wvl),3,20))
	#
	#
	if len(ind_i) >= 1:
			cof=polyfit(wvl[0:nel],specfl[0:nel],1)
			#
			#	   This makes the transition smooth, if not
			#	   there is a break as the linfit is evaluated
			#	   not only at the last point and the coordinate
			#	   at the origin of the fit makes the extrapolated
			#	   curve jump at this point
			#
			if hasattr(new_wvl, '__getitem__'):
				temp=polyval(cof,new_wvl[ind_i])
			else:
				temp=polyval(cof,new_wvl)

			result[ind_i]=(temp+(specfl[0]-temp[-1])).clip(0)
	#### FIXME ####
	if len(ind_f) >= 1:
		print ("DEBUG==============", ind_f)
		cof=polyfit(wvl[-1*nel:],specfl[-1*nel:],1)
		if not hasattr(new_wvl, '__getitem__'):
			new_wvl = [new_wvl]
		temp=polyval(cof,new_wvl[ind_f[0]])
		if not hasattr(temp, '__getitem__'):
			temp = [temp]
		result[ind_f[0]]=(temp+(specfl[-1]-temp[0])).clip(0)
		print (new_wvl, ind_f, temp, specfl)
	#### END FIXME ####
	print (result)
	if working_unit == None:
		return result
	else:
		return result * working_unit
	#

def e_to_ph(value):
	return value*u.ph/u.electron

def interpol2newwave(spec,wave,newwave):
 
   #return np.interp(wave,self.wave,self.response)
   same_unit = newwave.unit     
   _xnew = newwave.value
   _xx = wave.to(same_unit).value
   if hasattr(spec,'unit'):
       _yy = spec.value
       spec_unit = spec.unit     
   else:
       _yy = spec
       spec_unit = 1
         
   _ynew = np.interp(_xnew,_xx,_yy) 
   return _ynew*spec_unit 


def average_bandpass(spec,bandpass,wave):
 
   #return np.interp(wave,self.wave,self.response)
   same_unit = wave.unit     
   _xx = wave.value
   if hasattr(spec,'unit'):
       _yy = spec.value * bandpass
       spec_unit = spec.unit     
   else:
       _yy = spec * bandpass
       spec_unit = 1
         
   avg_spec = scipy.integrate.simps(_yy,_xx)/scipy.integrate.simps(bandpass,_xx)
            
   return avg_spec*spec_unit 
