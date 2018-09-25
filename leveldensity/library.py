################################################################################
# Author: Alec Golas                                                           #
# Date: September 25th, 2018                                                   #
# Title: library.py                                                            #
# Default parameters and load nuclei data                                      #
################################################################################

import os
import sys
import json

################################################################################
class Loader:
    """This class is intended to load the input specifications from a python
    dictionary into the proper variables in order to load the correct target data parameters to
    perform the level density calculations"""

    def __init__(self, inputdict):
        self.parameters = self.load(inputdict)



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
        return calc_parameters

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
############################# End of library.py ################################
