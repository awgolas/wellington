################################################################################
# Author: Alec Golas                                                           #
# Date: September 18th, 2018                                                   #
# Title: leveldensity.py                                                       #
# Determines cr52 Level Densities for plotting using different models          #
################################################################################
from __future__ import print_function
import matplotlib
matplotlib.use('agg')

import os
import sys
import numpy as np
import math
import json
import matplotlib.pyplot as plt
from utilities import Math
from library import Loader

import level_parameters

def run():

    excitation_energy = np.linspace(2,8,num=50)
    parity = 1
    target = '52Cr'
    temp = 8
    spin = [1.0]
    for j in spin:
        input = {'target' : target,
                 'spin'   : j,
                 'parity' : parity,
                 'excitation_energy' : excitation_energy}

    gcm = level_parameters.CompositeGilbertCameronParameters(input).mass
    bgsm = level_parameters.BackShiftedFermiGasParameters(input).leveldensity
    #print(bgsm)
    print(gcm)
    #p = plt.scatter(excitation_energy, bgsm)
    #p = plt.scatter(excitation_energy, gcm)
    #p.figure.savefig('levden_comp.png')
        #gcm = GilbertCameronModel(input)

        #rho_fgm = bsfgm.leveldensity()
        #rho_gcm = gcm.leveldensity()

# ################################################################################
# class BackShiftedFermiGasModel(level_parameters.BackShiftedFermiGasParameters):
#
#     def leveldensity(self):
#
#         j = self.spin
#         pi = self.pi
#         energy=self.excitation_energy
#
#         eff_energy = self.eff_energy
#         a = self.afgm
#
#
#         sigma2 = self.sigma_squared
#
#         rho0  = (0.5)
#         rho1  =       (2.0*j+1.0)
#         rho2  =       1.0/(2.0*(2.0*pi)**0.5*sigma2**(1.5))
#         rho3  =       np.exp(((j+0.5)**2.0)*(1/(2*sigma2)))
#         rho4  =       (pi**0.5)/12.0
#         rho5  =       np.exp(2.0*(a*eff_energy)**0.5)
#         rho6  =       1.0/(a**0.25*eff_energy**1.25)
#
#         rho = rho1*rho2*rho3*rho4*rho5*rho6
#         return rho
#
#     def cum_level_density(self, leveldensity):
#         util = Utilities()
#         cld = util.summation(leveldensity)
#         return cld
#
# ################################################################################
# class GilbertCameronModel(level_parameters.CompositeGilbertCameronParameters):
#     """ Utilizes the Composite Gilbert Cameron Model (GCM) in order to determine
#     level density of a compound nucleus as a function of energy, spin, and parity.
#     Subfunctions of GCM (e.g., cutoff_energy, temp, delta, and matching_energy)
#     are included in CompositeGilbertCameronParameters class in level_param.py"""
#
#
#     @property
#     def leveldensity(self):
#
#         cld = self.cum_level_density
#         temp = self.temperature
#         rho = cld/temp
#         return rho
#
#     @property
#     def cum_level_density(self):
#         temp = self.temperature
#         excitation_energy = self.excitation_energy
#         cutoff_energy = self.cutoff_energy
#
#         exponent = (excitation_energy - cutoff_energy)/temp
#         cld = np.exp(exponent)
#         return cld
#
#
# ################################################################################
""" This section should calculate the nuclear deformation values"""

"""dynamic deformation

    a2dyn = b*(-1.25*y/(1-x))

    where:
        y -> angular momentum parameter
        x -> fissility parameter
        b -> adjustable parameter

    y = 1.9249*I(I+1)*(I*(I+1)/(eta*A^(7/3)))

    where:
        I -> angular momentum
        eta -> neutron proton difference term
        A -> number of nucleons in nuclei

    eta = 1 - 1.7826*(N-Z)^2*A^-2

    ##

    a2(T,I) = a_{g.s.}*h(T) + a2dyn"""


################################################################################
if __name__ == "__main__":
    run()
################################################################################
######################## End of leveldensity.py ################################
