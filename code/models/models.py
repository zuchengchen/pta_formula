from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from enterprise.signals import parameter
from enterprise.signals import selections
from enterprise.signals import signal_base
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals
from enterprise.signals import utils

from utils import get_noise_dict

# import omegagw_funcs as OG
# from omegagw_funcs import Ogw2

H0 = 2.192711267238057e-18
from numpy import pi
from omegagw_funcs import STT, SGeneral, SInflation

# range for log10@A
logA_min, logA_max = -20, -10

@signal_base.function
def GWTT(f, log10_A, kappa):
    df = np.diff(np.concatenate((np.array([0]), f[::2])))
    return STT(f, 10**log10_A, kappa) * np.repeat(df, 2)

@signal_base.function
def GWGeneral(f, log10_A, kappa):
    df = np.diff(np.concatenate((np.array([0]), f[::2])))
    return SGeneral(f, 10**log10_A, kappa) * np.repeat(df, 2)

@signal_base.function
def GWInflation(f, log10_A):
    df = np.diff(np.concatenate((np.array([0]), f[::2])))
    return SInflation(f, 10**log10_A) * np.repeat(df, 2)


#### Extra model componentents not part of base enterprise ####

@signal_base.function
def free_spectrum(f, log10_rho=None):
    """
    Free spectral model. PSD  amplitude at each frequency
    is a free parameter. Model is parameterized by
    S(f_i) = \rho_i^2 * T,
    where \rho_i is the free parameter and T is the observation
    length.
    """
    return np.repeat(10**(2*log10_rho), 2)


#### Model component building blocks ####

def white_noise_block(vary=False):
    """
    Returns the white noise block of the model:

        1. EFAC per backend/receiver system
        2. EQUAD per backend/receiver system
        3. ECORR per backend/receiver system

    :param vary:
        If set to true we vary these parameters
        with uniform priors. Otherwise they are set to constants
        with values to be set later.
    """

    # define selection by observing backend
    selection = selections.Selection(selections.by_backend)

    # white noise parameters
    if vary:
        efac = parameter.Uniform(0.01, 10.0)
        equad = parameter.Uniform(-8.5, -5)
        ecorr = parameter.Uniform(-8.5, -5)
    else:
        efac = parameter.Constant()
        equad = parameter.Constant()
        ecorr = parameter.Constant()

    # white noise signals
    ef = white_signals.MeasurementNoise(efac=efac, selection=selection)
    eq = white_signals.EquadNoise(log10_equad=equad, selection=selection)
    ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=selection)

    # combine signals
    s = ef + eq + ec

    return s


def red_noise_block(prior='log-uniform', Tspan=None):
    """
    Returns red noise model:

        1. Red noise modeled as a power-law with 30 sampling frequencies

    :param prior:
        Prior on log10_A. Default if "log-uniform". Use "uniform" for
        upper limits.
    :param Tspan:
        Sets frequency sampling f_i = i / Tspan. Default will
        use overall time span for indivicual pulsar.

    """

    # red noise parameters
    if prior == 'uniform':
        log10_A = parameter.LinearExp(-20, -11)
    elif prior == 'log-uniform':
        log10_A = parameter.Uniform(-20, -11)
    else:
        raise ValueError('Unknown prior for red noise amplitude!')

    gamma = parameter.Uniform(0, 7)

    # red noise signal
    pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
    rn = gp_signals.FourierBasisGP(pl, components=50, Tspan=Tspan)

    return rn


def common_red_noise_block(psd='inflation', 
                           prior='log-uniform',
                           Tspan=None, 
                           orf="TT",
                           name=None,
                           kappa_val=None,
                           logA_min=-30, 
                           logA_max=-8):
    """
    Returns common red noise model:

        1. Red noise modeled with user defined PSD with
        30 sampling frequencies. Available PSDs are
        ['powerlaw', 'turnover' 'spectrum']

    :param psd:
        PSD to use for common red noise signal. Available options
        are ['powerlaw', 'turnover' 'spectrum']
    :param prior:
        Prior on log10_A. Default if "log-uniform". Use "uniform" for
        upper limits.
    :param Tspan:
        Sets frequency sampling f_i = i / Tspan. Default will
        use overall time span for indivicual pulsar.
    :param orf:
        String representing which overlap reduction function to use.
        By default we do not use any spatial correlations. Permitted
        values are ['hd', 'dipole', 'monopole'].
    :param name: Name of common red process

    """

    orfs = {'TT': utils.TT_orf(), 'ST': utils.ST_orf(),
            'VL': utils.VL_orf(), 'SL':utils.SL_orf()}

    if name==None: 
        name = orf
    
    Agw_name = 'log10_A_{}'.format(name)
    if prior == 'uniform':  
        log10_Agw = parameter.LinearExp(logA_min, logA_max)(Agw_name)
    elif prior == 'log-uniform':
        log10_Agw = parameter.Uniform(logA_min, logA_max)(Agw_name)   
        
        
    if psd in ['SMBH-TT', 'SMBH-General']:
        kappa_name = 'kappa_{}'.format(name)
        kappa_gw = parameter.Uniform(0, 10)(kappa_name)
            
        if psd == "SMBH-TT":
            cpl = GWTT(log10_A=log10_Agw, kappa=kappa_gw)
        if psd == "SMBH-General":
            cpl = GWGeneral(log10_A=log10_Agw, kappa=kappa_gw)      
         
        
    if psd == 'inflation':
        cpl = GWInflation(log10_A=log10_Agw)

    crn = gp_signals.FourierBasisCommonGP(cpl, orfs[orf], components=30, Tspan=Tspan, name=name)

    return crn

def modelInflation(psrs, upper_limit=False, bayesephem=True):
    """
    Reads in list of enterprise Pulsar instance and returns a PTA
    instantiated with model 3A from the analysis paper:

    per pulsar:
        1. fixed EFAC per backend/receiver system
        2. fixed EQUAD per backend/receiver system
        3. fixed ECORR per backend/receiver system
        4. Red noise modeled as a power-law with 30 sampling frequencies
        5. Linear timing model.

    global:
        1. GWB with HD correlations modeled with user defined PSD with
        30 sampling frequencies. Available PSDs are
        ['powerlaw', 'turnover' 'spectrum']
        2. Optional physical ephemeris modeling.

    :param psd:
        PSD to use for common red noise signal. Available options
        are ['powerlaw', 'turnover' 'spectrum'] 'powerlaw' is default
        value.
    :param gamma_common:
        Fixed common red process spectral index value. By default we
        vary the spectral index over the range [0, 7].
    :param upper_limit:
        Perform upper limit on common red noise amplitude. By default
        this is set to False. Note that when perfoming upper limits it
        is recommended that the spectral index also be fixed to a specific
        value.
    :param bayesephem:
        Include BayesEphem model. Set to False by default
    """

    amp_prior = 'uniform' if upper_limit else 'log-uniform'

    # find the maximum time span to set GW frequency sampling
    tmin = [p.toas.min() for p in psrs]
    tmax = [p.toas.max() for p in psrs]
    Tspan = np.max(tmax) - np.min(tmin)

    # white noise
    s = white_noise_block(vary=False)

    # red noise
    s += red_noise_block(prior=amp_prior, Tspan=Tspan)

    # common red noise block
    s += common_red_noise_block(psd='inflation', prior=amp_prior, Tspan=Tspan, orf='TT')
    s += common_red_noise_block(psd='inflation', prior=amp_prior, Tspan=Tspan, orf='ST')    
    s += common_red_noise_block(psd='inflation', prior=amp_prior, Tspan=Tspan, orf='VL')
    s += common_red_noise_block(psd='inflation', prior=amp_prior, Tspan=Tspan, orf='SL')

    # ephemeris model
    if bayesephem:
        s += deterministic_signals.PhysicalEphemerisSignal(use_epoch_toas=True)

    # timing model
    s += gp_signals.TimingModel()

    # set up PTA
    pta = signal_base.PTA([s(psr) for psr in psrs])

    # set white noise parameters
    noisedict = get_noise_dict(psrlist=[p.name for p in psrs])
    pta.set_default_params(noisedict)

    return pta

def modelSMBH(psrs, upper_limit=False, bayesephem=True):
    """
    Reads in list of enterprise Pulsar instance and returns a PTA
    instantiated with model 3A from the analysis paper:

    per pulsar:
        1. fixed EFAC per backend/receiver system
        2. fixed EQUAD per backend/receiver system
        3. fixed ECORR per backend/receiver system
        4. Red noise modeled as a power-law with 30 sampling frequencies
        5. Linear timing model.

    global:
        1. GWB with HD correlations modeled with user defined PSD with
        30 sampling frequencies. Available PSDs are
        ['powerlaw', 'turnover' 'spectrum']
        2. Optional physical ephemeris modeling.

    :param psd:
        PSD to use for common red noise signal. Available options
        are ['powerlaw', 'turnover' 'spectrum'] 'powerlaw' is default
        value.
    :param gamma_common:
        Fixed common red process spectral index value. By default we
        vary the spectral index over the range [0, 7].
    :param upper_limit:
        Perform upper limit on common red noise amplitude. By default
        this is set to False. Note that when perfoming upper limits it
        is recommended that the spectral index also be fixed to a specific
        value.
    :param bayesephem:
        Include BayesEphem model. Set to False by default
    """

    amp_prior = 'uniform' if upper_limit else 'log-uniform'

    # find the maximum time span to set GW frequency sampling
    tmin = [p.toas.min() for p in psrs]
    tmax = [p.toas.max() for p in psrs]
    Tspan = np.max(tmax) - np.min(tmin)

    # white noise
    s = white_noise_block(vary=False)

    # red noise
    s += red_noise_block(prior=amp_prior, Tspan=Tspan)

    # common red noise block
    s += common_red_noise_block(psd='SMBH-TT', prior=amp_prior, Tspan=Tspan, orf='TT', logA_max=-13)  
    s += common_red_noise_block(psd='SMBH-General', prior=amp_prior, Tspan=Tspan, orf='ST', logA_max=-13)  
    s += common_red_noise_block(psd='SMBH-General', prior=amp_prior, Tspan=Tspan, orf='VL', logA_max=-13)  
    s += common_red_noise_block(psd='SMBH-General', prior=amp_prior, Tspan=Tspan, orf='SL', logA_max=-13)  

    # ephemeris model
    if bayesephem:
        s += deterministic_signals.PhysicalEphemerisSignal(use_epoch_toas=True)

    # timing model
    s += gp_signals.TimingModel()

    # set up PTA
    pta = signal_base.PTA([s(psr) for psr in psrs])

    # set white noise parameters
    noisedict = get_noise_dict(psrlist=[p.name for p in psrs])
    pta.set_default_params(noisedict)

    return pta