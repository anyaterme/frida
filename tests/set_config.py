import scipy.constants as pconst

def get_filter(filter='H'):
    #
    #       Gets the characteristics of filter: cut and returns it as a 
    #       speccurve object
    #
    if (filter == 'H'):
       name='H'
       cut_on=1.490  # microns
       cut_off=1.780
       trans=0.85
    elif (filter == 'J'):
       name='J'
       cut_on=1.170
       cut_off=1.330
       trans=0.85
    elif (filter == 'Ks'):
       name='Ks'
       cut_on=1.990
       cut_off=2.310
       trans=0.85
    elif (filter == 'Jc'):
       name='Jc'
       cut_on=1.216
       cut_off=1.244
       trans=0.8

    lambda_eff=(cut_on+cut_off)/2.
    bandwidth=(cut_off-cut_on)
    energy_phot = pconst.h*pconst.c/(lambda_eff*1.e-6) ## in SI units

    #
    return {'name':name,'lamb_eff':lambda_eff,'bandwidth':bandwidth,'trans':trans,'energy_phot':energy_phot}

def par_gtcao(filter='H'):

    if (filter == 'H'):
       name='H'
       throughput=0.85
    elif (filter == 'J'):
       name='J'
       throughput=0.85
    elif (filter == 'Ks'):
       name='Ks'
       troughput=0.85
    elif (filter == 'Jc'):
       name='Jc'
       throughput=0.85

    return {'name':name,'throughput':throughput}


def par_instru(filter='H'):
        
    ron = 15.   ## in e-
    gain = 3.5
    welldepth = 2.e5  ## in e-
    dark = 0.1 ## dark current [e-/sec]

    if (filter == 'H'):
       name='H'
       qe=0.9
       effic_cam=0.85
       effic_relay=0.85
    elif (filter == 'J'):
       name='J'
       qe=0.9
       effic_cam=0.85
       effic_relay=0.85
    elif (filter == 'Ks'):
       name='Ks'
       qe=0.9
       effic_cam=0.85
       effic_relay=0.85
    elif (filter == 'Jc'):
       name='Jc'
       qe=0.9
       effic_cam=0.85
       effic_relay=0.85
    #
    return {'name':name,'qe':qe,'effic_cam':effic_cam,'effic_relay':effic_relay,'ron':ron,'welldepth':welldepth,'dark':dark}

def par_atm(filter='H'):
        
    if (filter == 'H'):
       name='H'
       trans=1.0
       mag_sqarc=13.6
    elif (filter == 'J'):
       name='J'
       trans=1.0
       mag_sqarc=17.
    elif (filter == 'Ks'):
       name='Ks'
       trans=0.9
       mag_sqarc=13.
    elif (filter == 'Jc'):
       name='Jc'
       trans=1.0
       mag_sqarc=17.
    #
    return {'name':name,'trans':trans,'mag_sqarc':mag_sqarc}

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
    
