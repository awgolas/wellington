################################################################################
# Author: Alec Golas                                                           #
# Date: September 24th, 2018                                                   #
# Title: level_param.py                                                        #
# Determines level density subfunction values                                  #
################################################################################
from __future__ import print_function
import matplotlib

import os
import sys
import numpy as np
import math
from utilities import Utilities
from leveldensity import NuclearData

################################################################################
class GeneralParameters(NuclearData):
    @property
    def shell_correction(self):

        mass = self.nuc_mass
        A    = self.mass_number
        N    = self.num_neutrons
        Z    = self.num_protons

        k  = 1.79
        a  = [15.677, 18.56]
        ci = (1-k*((N-Z)/A)**2)
        c  = [a[0]*ci, a[1]*ci, 0.717, 1.21129]

        evol  = -c[0]*A
        esur  = c[1]*A**(0.66667)
        ecoul = c[2]*(Z**2.0/(A**(0.33333))) - c[3]*(Z**2.0/A)

        if A%2 == 0:
            if Z%2 == 0:
                delta_m = -11.0/(A**(0.5))
            else:
                delta_m = 11.0/(A**(0.5))
        else:
            delta_m = 0.0

        m_n = 8.07144
        m_h = 7.28899

        m_ldm = m_n*N + m_h*Z + evol + esur + ecoul + delta_m

        delta_w = mass - m_ldm

        return delta_w

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

        a_param = atilda*(1.0 + u*(1.0/delta_w)*(1.0 - self.exp**(-1*gamma*u)))
        return a_param

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
        a   = atilda*(1.0 + delta_w/u*(1.0 - np.exp(-1.0*gamma*u)))
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

        s_squared[lindex] = sigma_d2 + excitation_energy[lindex]*(1.0/s_n)\
            *(spin_cutoff_bn - sigma_d2)
        s_squared[uindex] = self.sigma_f2(excitation_energy[uindex]))

        return s_squared

    def sigma_f2(self, excitation_energy):

        A       = self.mass_number
        delta_w = self.shell_correction
        gamma   = self.gamma
        atilda  = self.atilda
        delta   = self.delta

        u   = excitation_energy - delta
        a   = atilda*(1.0 + delta_w/u*(1.0 - self.exp**(-1.0*gamma*u)))
        sf2 = 0.01389*A**(1.66667)/atilda*((a*u)**(0.5))
        return sf2




################################################################################
class CompositeGilbertCameronParameters(GeneralParameters):


    @property
    def log_rho(self):
        excitation_energy = self.excitation_energy
        fgm_rho = fgm.leveldensity(energy=excitation_energy)
        ln_rho = np.log(fgm_rho)

        return ln_rho

    @property
    def cutoff_energy(self):
        util = Utilities()

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
        excitation_energy = self.excitation_energy

        e_m = 2.33 + 253/A + delta
        return e_m

    @property
    def temperature(self):
        util = Utilities()
        e_m = self.matching_energy
        excitation_energy = self.excitation_energy
        ln_rho = self.log_rho

        dfdx_ln_rho = util.dfdx_1d(excitation_energy, ln_rho)

        em_index = util.nearest_value_index(excitation_energy, e_m)

        temp = 1.0/dfdx_ln_rho[em_index]
        return temp
################################################################################
