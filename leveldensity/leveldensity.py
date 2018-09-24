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
from utilities import Utilities

from level_param import CompositeGilbertCameronParameters, BackShiftedFermiGasParameters

def run():

    excitation_energy = list(np.linspace(2,10,num=100))
    parity = 1
    target = '51Cr'
    temp = 8
    spin = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    for j in spin:
        input = {'target' : target,
                 'spin'   : j,
                 'parity' : parity,
                 'excitation_energy' : excitation_energy}


        rho = BackShiftedFermiGasModel().run(input)
        label = 'Spin:{}'.format(j)
        p = plt.plot(excitation_energy, rho, label=label)
    plt.ylabel('Level Density (1/MeV)')
    plt.xlabel('Excitation Energy (MeV)')
    plt.legend(loc=0)
    plt.title('LD Calculation of 51Cr using Back Shifted Fermi Gas Model')
    plt.savefig('levden51.png')

    #BackShiftedFermiGasModel().main()
################################################################################
class Loader:
    """This class is intended to load the input specifications from a python
    dictionary into the proper variables in order to load the correct target data parameters to
    perform the level density calculations"""



    def load(self, inputdict):
        input_params = self.input(inputdict)
        target = input_params['target']
        nuclear_data = self.library(target)
        calc_parameters = {'target' : target,
                           'spin'   : input_params['spin'],
                           'parity' : input_params['parity'],
                           'excitation_energy' : input_params['excitation_energy'],
                           'A'               : nuclear_data['A'],
                           'Z' : nuclear_data['Z'],
                           'mass' : nuclear_data['mass'],
                           'Bn' : nuclear_data['Bn']}
        with open('parameters.json', 'w') as f:
            json.dump(calc_parameters, f)

    def input(self, inputdict):
        output = {}
        default = {'target' : None,
                   'spin'   : 0,
                   'parity' : 0,
                   'temperature' : 8,
                   'excitation_energy' : 0}

        for parameter in default.keys():
            try:
                val = inputdict[parameter]
            except:
                val = default[parameter]

            output[parameter] = val

        return output


    def library(self, target):
        with open('data.json', 'r') as f:
            data = json.load(f)
        isotope_data = data[target]
        return isotope_data
################################################################################
class NuclearData():
    """Collects the formatted input data and performs some simple operations
    to be used in the level density calculations"""

    def parameters(self, arg):
        with open('parameters.json', 'r') as f:
            data = json.load(f)
        out = data[arg]
        return out

    @property
    def exp(self):
        return 2.718281828459

    @property
    def hbar(self):
        return 6.58211928e-16

    @property
    def spin(self):
        j = self.parameters('spin')
        return j

    @property
    def pi(self):
        parity = self.parameters('parity')
        return parity

    @property
    def excitation_energy(self):
        ex_energy = np.asarray(self.parameters('excitation_energy'))
        return ex_energy

    @property
    def num_protons(self):
        num_p = self.parameters('Z')
        return num_p

    @property
    def num_neutrons(self):
        Z = self.parameters('Z')
        A = self.parameters('A')
        nn = A-Z
        return nn

    @property
    def mass_number(self):
        A = self.parameters('A')
        return A
    @property
    def nuc_mass(self):
        mass = self.parameters('mass')
        return mass

    @property
    def separation_energy(self):
        bn = self.parameters('Bn')
        return bn


################################################################################
class BackShiftedFermiGasModel(NuclearData):


    def leveldensity(self, energy=NuclearData().excitation_energy):

        j = self.spin
        pi = self.pi

        eff_energy = self.eff_energy
        a = self.afgm


        sigma2 = self.sigma_squared

        rho0  = (0.5)
        rho1  =       (2.0*j+1.0)
        rho2  =       1.0/(2.0*(2.0*pi)**0.5*sigma2**(1.5))
        rho3  =       self.exp**(((j+0.5)**2.0)*(1/(2*sigma2)))
        rho4  =       (pi**0.5)/12.0
        rho5  =       self.exp**(2.0*(a*eff_energy)**0.5)
        rho6  =       1.0/(a**0.25*eff_energy**1.25)

        rho = rho1*rho2*rho3*rho4*rho5*rho6
        return rho

    def cum_level_density(self, leveldensity):
        util = Utilities()
        cld = util.summation(leveldensity)
        return cld




################################################################################
class GilbertCameronModel(CompositeGilbertCameronParameters):
    """ Utilizes the Composite Gilbert Cameron Model (GCM) in order to determine
    level density of a compound nucleus as a function of energy, spin, and parity.
    Subfunctions of GCM (e.g., cutoff_energy, temp, delta, and matching_energy)
    are included in CompositeGilbertCameronParameters class in level_param.py"""

    @property
    def leveldensity(self):

        cld = self.cum_level_density
        temp = self.temperature
        rho = cld/temp
        return rho

    @property
    def cum_level_density(self):
        temp = self.temperature
        excitation_energy = self.excitation_energy
        cutoff_energy = self.cutoff_energy

        exponent = (excitation_energy - cutoff_energy)/temp
        cld = np.exp(exponent)
        return cld


################################################################################
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
