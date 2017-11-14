import numpy as np
from math import ceil

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
    pulse_width_units = ceil(pulse_width / delt_wave)
    pulsodelta = np.zeros(3*pulse_width_units)
    for i in range(pulse_width_units):
        pulsodelta[i+pulse_width_units] = 1
    pulsodelta = pulsodelta / pulsodelta.sum()
    conv_response = np.convolve(response,pulsodelta,'same')
    return conv_response

def interpolate(wvl,specfl,new_wvl,unity='None'):
    #
    #       Interpolation to a given wvl array. Returned values
    #       are clipped to 0, and if it is a transmission curve,
    #       also to ones.
    #
    #       Order=1 is used, as larger orders produce
    #       weird effects around corners.
    #
    from scipy.interpolate import InterpolatedUnivariateSpline as inter
    from scipy import polyval,polyfit
    #
    inter_func=inter(wvl,specfl,k=1)
    if (unity == 'perone'):
            result=inter_func(new_wvl).clip(0,1)
    elif (unity == 'percent'):
            result=inter_func(new_wvl).clip(0,100)
    else:
            result=inter_func(new_wvl)
    #
    #       To avoid weird extrapolations, a lineal fit to the
    #       extreme points of the original curve is used
    #
    ind_i=np.where(new_wvl <= wvl.min())[0]
    ind_f=np.where(new_wvl >= wvl.max())[0]
    nel=np.clip(0.1*len(wvl),3,20)
    #
    #
    if len(ind_i) >= 1:
            cof=polyfit(wvl[0:nel],specfl[0:nel],1)
            #
            #       This makes the transition smooth, if not
            #       there is a break as the linfit is evaluated
            #       not only at the last point and the coordinate
            #       at the origin of the fit makes the extrapolated
            #       curve jump at this point
            #
            temp=polyval(cof,new_wvl[ind_i])
            result[ind_i]=(temp+(specfl[0]-temp[-1])).clip(0)
    #
    if len(ind_f) >= 1:
            cof=polyfit(wvl[-1*nel:],specfl[-1*nel:],1)
            temp=polyval(cof,new_wvl[ind_f])
            temp2=polyval(cof,new_wvl[-1])
            result[ind_f]=(temp+(specfl[-1]-temp[0])).clip(0)
    #
    return result
    #
