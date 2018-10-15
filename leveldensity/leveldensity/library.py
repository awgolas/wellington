################################################################################
# Author: Alec Golas                                                           #
# Date: September 25th, 2018                                                   #
# Title: library.py                                                            #
# Default parameters and load nuclei data                                      #
################################################################################

import os
import sys
import json
import re
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

        calc_parameters = {'target_label'      : target,
                           'compound_label'    : nuclear_data['compound'],
                           'parity'            : input_params['parity'],
                           'spin'              : nuclear_data['spin'],
                           'excitation_energy' : input_params['excitation_energy'],
                           'A'                 : nuclear_data['A'],
                           'Z'                 : nuclear_data['Z'],
                           'mass'              : nuclear_data['mass'],
                           'shell_correction'  : nuclear_data['delta_w'],
                           'Bn'                : nuclear_data['Io'],
                           'D0'                : nuclear_data['D0'],
                           'D0_err'            : nuclear_data['D0_err']}
        return calc_parameters

    def input(self, inputdict):
        output = {}
        default = {'target' : None,
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

        r = re.compile("([0-9]+)([a-zA-Z]+)")
        m = r.match(target)
        compound_z = str(int(m.group(1)) + 1)
        compound_nucleus = compound_z + m.group(2)

        target_data = data[target]
        compound_data = data[compound_nucleus]

        nuclear_data = {'target'   : target,
                        'compound' : compound_nucleus,
                        'A'        : compound_data['A'],
                        'Z'        : compound_data['Z'],
                        'mass'     : compound_data['mass'],
                        'delta_w'  : compound_data['shell_correction'],
                        'spin'     : target_data['Io'],
                        'Bn'       : target_data['Bn'],
                        'D0'       : target_data['D0'],
                        'D0_err'   : target_data['dD0']}
        return nuclear_data

################################################################################
class Parameters(Loader):

    """
    Instantiates the dictionary of problem parameters into class methods which
    is faster and simpler than instantiating the Loader class each time a
    parameter is required.
    """

    def __init__(self, inputdict):
        Loader.__init__(self, inputdict)


    @property
    def target_label(self):
        return self.parameters['target_label']

    @property
    def compound_label(self):
        return self.parameters['compound_label']

    @property
    def pi(self):
        return self.parameters['parity']

    @property
    def spin(self):
        return self.parameters['spin']

    @property
    def excitation_energy(self):
        return self.parameters['excitation_energy']

    @property
    def mass_number(self):
        return self.parameters['A']

    @property
    def num_protons(self):
        return self.parameters['Z']

    @property
    def mass(self):
        return self.parameters['mass']

    @property
    def shell_correction(self):
        return self.parameters['shell_correction']

    @property
    def separation_energy(self):
        return self.parameters['Bn']

    @property
    def res_separation_spacing(self):
        return self.parameters['D0']

    @property
    def res_separation_spacing_error(self):
        return self.parameters['D0_err']


################################################################################    
class Output:
    pass
################################################################################
############################# End of library.py ################################
