################################################################################
# Author: Alec Golas                                                           #
# Date: September 24th, 2018                                                   #
# Title: level_parameters.py                                                   #
# Determines level density subfunction values                                  #
################################################################################
from __future__ import print_function

import os
import sys
import numpy as np
import math
from utilities import Math
from library import Parameters

################################################################################
class GeneralParameters(Parameters):

    """ 
    Calculates the general level density parameters of the nuclei being observed
    that are general to all level density models.
    """

    @property
    def num_neutrons(self):
        Z = self.num_protons
        A = self.mass_number
        return A - Z

    @property
    def delta(self):
        A = self.mass_number
        Z = self.num_protons

        if A%2 == 0:
            if Z%2 == 0:
                n = 2
            else:
                n = 0
        else:
            n = 1

        d = 12/(A**(0.5))
        return d

    @property
    def global_spin_cutoff(self):
        A = self.mass_number

        sigma_d2 = (0.83*A**0.26)**2.0
        return sigma_d2

################################################################################
class BackShiftedFermiGasParameters(GeneralParameters):

    """
    Level density model that utilizes Fermi-Gas Theory as defined in the RIPL
    Reference Input Parameter Library.

    Reference:
    Capote, R., et al. "RIPL - Reference Input Parameter Library for Calculation
        of Nuclear Reactions and Nuclear Data Evaluations."
        Nuclear Data Sheets, vol. 110, no. 12, 2009, pp. 3107 - 3214.,
        doi:10.1016/j.nds.2009.10.004.
    """


    @property
    def eff_energy(self):
        excitation = self.excitation_energy
        pairing    = self.fgm_delta

        u = excitation - pairing
        return u

    @property
    def fgm_delta(self):
        delta      = self.delta
        correction = 0.173015

        d = delta + correction
        return d

    @property
    def afgm(self):
        atilda  = self.atilda
        delta_w = self.shell_correction
        u       = self.eff_energy
        gamma   = self.gamma

        f_u = 1 - np.exp(-1.0*gamma*u)

        a_param = atilda*(1.0 + delta_w/u*f_u)
        return a_param


    @property
    def atilda(self):
        A = self.mass_number
        alpha = 0.0722396
        beta = 0.195267

        at = alpha*A + beta*A**(0.66667)

        return at

    @property
    def gamma(self):
        gam0 = 0.410289
        A    = self.mass_number

        g = gam0/(A**(0.33333))

        return g

    @property
    def spin_cutoff_bn(self):
        ex_engy = self.separation_energy
        A       = self.mass_number
        delta_w = self.shell_correction
        gamma   = self.gamma
        atilda  = self.atilda
        delta   = self.fgm_delta

        u   = ex_engy - delta
        f_u = 1.0 - np.exp(-1.0*gamma*u)
        a = atilda*(1.0 + delta_w/u*f_u)
        sf2 = 0.01389*A**(1.66667)/atilda*((a*u)**(0.5))
        return sf2

    @property
    def sigma_squared(self):
        s_n = self.separation_energy
        excitation_energy = self.excitation_energy
        spin_cutoff_bn = self.spin_cutoff_bn
        sigma_d2 = self.global_spin_cutoff

        s_squared = np.zeros(np.shape(excitation_energy))

        uindex = np.where(excitation_energy >= s_n)
        lindex = np.where(excitation_energy < s_n)

        print(excitation_energy)
        s_squared[lindex] = sigma_d2 + excitation_energy[lindex]*(1.0/s_n)\
                *(spin_cutoff_bn - sigma_d2)
        s_squared[uindex] = self.sigma_f2(excitation_energy[uindex])

        return s_squared

    def sigma_f2(self, excitation_energy):

        A       = self.mass_number
        delta_w = self.shell_correction
        gamma   = self.gamma
        atilda  = self.atilda
        delta   = self.delta

        u   = excitation_energy - delta
        f_u = 1 - np.exp(-1.0*gamma*u)
        a = atilda*(1.0 + delta_w/u*f_u)
        sf2 = 0.01389*A**(1.66667)/atilda*((a*u)**(0.5))
        return sf2

    @property
    def rho_energy(self):

        a = self.afgm
        u = self.eff_energy
        pi = self.pi

        rho_e = np.exp(2.0*(a*u)**0.5)*pi**(0.5)/(12.0*a**(0.25)*u**(1.25))
        return rho_e

    @property
    def rho_jpi(self):
        j = self.spin
        pi = self.pi
        sigma2 = self.sigma_squared

        rho_jp = np.exp(-1.0*(j+0.5)**(2.0)/(2.0*sigma2))\
            *(0.5)*(2.0*j+1)/((8.0*pi*sigma2**(1.5))**(0.5))
        return rho_jp

    @property
    def leveldensity(self):

        rho_e = self.rho_energy
        rho_jpi = self.rho_jpi

        rho = rho_e*rho_jpi
        return rho


################################################################################
class CompositeGilbertCameronParameters(GeneralParameters):

    """
    Level density model that utilizes the Back Shifted Fermi-Gas Model for a
    high-energy region and Constant-Temperature Model for a low-energy region.
    Uses parameterized values defined in the RIPL Reference Input Parameter
    Library. Same reference as BSFGM
    """

    def __init__(self, inputval):
        GeneralParameters.__init__(self, inputval)
        bsfgm = BackShiftedFermiGasParameters(inputval)
        self.fgm_rho_energy = bsfgm.rho_energy
        self.rho_jpi = bsfgm.rho_jpi

    @property
    def eff_energy(self):
        excitation = self.excitation_energy
        pairing    = self.delta

        u = excitation - pairing
        return u


    @property
    def log_rho(self):
        fgm_rho = self.fgm_rho_energy
        ln_rho = np.log(fgm_rho)

        return ln_rho

    @property
    def cutoff_energy(self):
        util = Math()
        excitation_energy = self.excitation_energy
        e_m = self.matching_energy
        temp = self.temperature
        ln_rho = self.log_rho

        em_index = util.nearest_value_index(excitation_energy, e_m)

        e0 = e_m - temp*(np.log(temp) + ln_rho[em_index])
        return e0

    @property
    def matching_energy(self):

        A = self.mass_number
        delta = self.delta

        e_m = 2.33 + 253.0/A + delta
        return e_m

    @property
    def temperature(self):
        util = Math()
        e_m = self.matching_energy
        excitation_energy = self.excitation_energy
        ln_rho = self.log_rho

        dfdx_ln_rho = util.dfdx_1d(excitation_energy, ln_rho)

        em_index = util.nearest_value_index(excitation_energy, e_m)

        temp = 1.0/dfdx_ln_rho[em_index]
        return temp

    @property
    def ct_rho_energy(self):

        cld = self.cum_rho_energy
        temp = self.temperature
        rho = cld/temp
        return rho

    @property
    def rho_energy(self):
        util = Math()
        excitation_energy = self.excitation_energy
        e_m = self.matching_energy
        ct_rho_energy = self.ct_rho_energy
        fgm_rho_energy = self.fgm_rho_energy

        em_index = util.nearest_value_index(excitation_energy, e_m)

        rho = np.zeros(np.shape(excitation_energy))

        rho[0:em_index] = ct_rho_energy[0:em_index]
        rho[em_index:] = fgm_rho_energy[em_index:]

        return rho

    @property
    def cum_rho_energy(self):
        temp = self.temperature
        excitation_energy = self.excitation_energy
        cutoff_energy = self.cutoff_energy

        exponent = (excitation_energy - cutoff_energy)/temp
        cld = np.exp(exponent)
        return cld

    @property
    def leveldensity(self):
        rho_e = self.rho_energy
        rho_jpi = self.rho_jpi

        rho = rho_e*rho_jpi
        return rho

################################################################################
class EmpireGilbertCameronParameters(GeneralParameters):

    """
    Uses Constant Temperature and Fermi Gas Models to calculate level
    densities but with the EMPIRE-derived parameter values rather than RIPL-
    derived parameters. Parameters are calculated in
    GilbertCameronFunctionParameters class.

    Reference:
    Herman, M., et al. EMPIRE-3.2 Malta Modular System for Nuclear Reaction
        Calculations and Nuclear Data Evaluation User's Manual. NNDC, Brookhaven
        National Laboratory, Upton, USA, 5 Aug. 2013.
    """

    def __init__(self, inputval, model):
        GeneralParameters.__init__(self, inputval)
        fp = GilbertCameronFunctionParameters(model)
        self.atilda_a = fp.atilde_a
        self.atilda_b = fp.atilda_b
        self.gamma = fp.gamma

        cgcp = CompositeGilbertCameronParameters(inputval)
        self.matching_energy = cgcp.matching_energy
        self.eff_energy = cgcp.eff_energy
        self.temperature = cgcp.temperature

    @property
    def atilda(self):
        A = self.mass_number
        a = self.atilda_a
        b = self.atilda_b

        return a*A + b[0]*A**b[1]

    @property
    def sigma_squared(self):
        A = self.mass_number
        u = self.eff_energy
        a = self.a

        sigma2 = 0.146*A**(0.6667)*(a*u)**0.5
        return sigma2

    @property
    def a(self):
        atilda  = self.atilda
        delta_w = self.shell_correction
        u       = self.eff_energy
        gamma   = self.gamma

        f_u = 1 - np.exp(-1.0*gamma*u)

        a_param = atilda*(1.0 + delta_w/u*f_u)
        return a_param


################################################################################
class GilbertCameronFunctionParameters:

    """
    Calculates parameters to be used in Empire Gilbert Cameron level density
    model. Uses Ignatyuk, Arthur, or Iljinov models as defined in EMPIRE manual.
    Same reference as EGCP.
    """

    def __init__(self, model):

        if model == 'ignatyuk':
            self._param = (0.154, (6.3e-5, 2), -0.054)
        elif model == 'arthur':
            self._param = (0.1375, (-8.36e-5, 2), -0.054)
        elif model == 'iljinov':
            self._param = (0.144, (9.8e-2, 2.0/3.0), -0.051)

    @property
    def atilde_a(self):
        return self._param[0]

    @property
    def atilde_b(self):
        return self._param[1]

    @property
    def gamma(self):
        return self._param[2]

# ################################################################################
# class GeneralizedSuperfluidParameters(GeneralParameters):
#
#     @property
#     def a_stndrd(self):
#
#     @property
#     def a_crit(self):
#
#     @property
#     def a_cum(self):
#
#     @property
#     def temp_crit(self):
#
#     @property
#     def temp_stndrd(self):
#
#     @property
#     def temp_cum(self):
#
#     @property
#     def det_crit(self):
#
#     @property
#     def det_stndrd(self):
#
#     @property
#     def
#
#
# ######################################
#     @property
#     def crit_temp(self):
#         pass
#
#     @property
#     def crit_eng(self):
#         pass
#
#     @property
#     def eff_energy(self):
#         pass
#
#     @property
#     def eng_shift(self):
#         pass
#
#     @property
#     def cond_eng(self):
#
#
# ################################################################################
# ##################### End of level_parameters.py ###############################
