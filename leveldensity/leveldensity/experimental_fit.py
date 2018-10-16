################################################################################
# Author: Alec Golas                                                           #
# Date: October 4th, 2018                                                      #
# Title: experimental_fit.py                                                   #
# Fits experimental level densities to curves                                  #
################################################################################

import os
import sys
import re

import numpy as np
from scipy import stats

from utilities import Math
from library import Parameters

################################################################################
class PostLevelDensityLoad(Parameters):
    """
    Pulls the levels for a nucleus from the RIPL or OSLO gamma emissions
    libraries and sorts them into CLD(E, J, pi) and CLD(E) data sets
    """

    def __init__(self, inputdict):
        Parameters.__init__(self, inputdict)
        self.exp_cld = self.get_data(self)
        

    @property
    def RIPL_dir(self):
        return os.path.join('/home/agolas', 'empire', 'RIPL', 'levels')

    def get_data(self, source='RIPL'):
        source='RIPL'
        nucleus = (self.compound_label, self.mass_number, self.num_protons)
        parity = self.pi
        spin  = self.spin

        if source == 'RIPL':
            data = self.ripl_data
        elif source == 'Oslo':
            data = self.oslo_data

        return data


    @property
    def ripl_data(self):
            
        leveldir = self.RIPL_dir
        
        label = self.compound_label
        r = re.compile("([0-9]+)([a-zA-Z]+)")
        m = r.match(label)
        end_z = str(int(m.group(1)) + 1)
        end_label = end_z + m.group(2)
        

        mass_number = self.mass_number
        charge_number = self.num_protons

        file_name = 'z{:03d}.dat'.format(charge_number)
        with open(leveldir +'/'+ file_name, 'r') as f:
            data_file = list(f)

        target = False
        for n, line in enumerate(data_file):
            
            if label in line:
                target = True
                init = n+1
            
            if target:
                if end_label in line:
                    end = n-1
                    break

        raw_nucleus_data = data_file[init:end]
        formatted_nucleus_data = []

        for line in raw_nucleus_data:
            if line.startswith(' '*31): 
                raw_nucleus_data.remove(line)
                continue

            clean_line = line.strip()
            formatted_nucleus_data.append(clean_line)

        jpi_col = []
        full_list = []
        cld_hist = {}
        for line in formatted_nucleus_data:
            sline = line.split()
            
            index = sline[0]
            eff_energy = float(sline[1])
            J = float(sline[2])
            pi = float(sline[3])

            Jpi = (J,pi)
            diag = (index, eff_energy, J, pi)
            full_list.append(diag)

            if Jpi not in jpi_col:
                jpi_col.append(Jpi)
                cld_hist[Jpi] = []

            cld_hist[Jpi].append(eff_energy)

        return cld_hist

    @property
    def Jpi(self):

        cld_hist = self.get_data()

        j_pis = list(cld_hist.keys())
        return j_pis

    def oslo_data(self, nucleus):
        #Do check to make sure nuclei is in the oslo folder
        #Open file with {nucleisymbol}.dat
        #That's basically it, they're preformatted
        #Also calculate CLD if it is requested
        #output data array as [level index, energy, level density]
        return 'Error my guy'

    def ripl_sort(self, **kwargs):

        #possibly accepted keywords:
            #spin
            #parity
        #calculates cumulative level density
        #calculates rho(E) from the derivative of the cumulative level density
        #outputs array as [excitation energy, spin, parity, rho]
        pass

################################################################################
class LevelDensityAnalyzer(PostLevelDensityLoad):

    def __init__(self, inputdict, JPi, calc_levels):
        PostLevelDensityLoad.__init__(self, inputdict)
        pldl = PostLevelDensityLoad(inputdict)
        self.exp_cld = pldl.exp_cld[JPi]
        self.JPi = JPi
        self.calc_levels = calc_levels

    @property
    def cld_hist(self):

        cld_energy = np.asarray(self.exp_cld)
        num_levels = len(cld_energy)
        index = np.arange(1, num_levels+1, step=1)

        histogram = np.array([index, cld_energy])
        return histogram

    @property
    def cld_equation_parameters(self):

        cld_linear = self.cld_hist
        j, pi = self.JPi


        cld_log = np.log(cld_linear)

        x,y = cld_log
        cld_loglinreg = stats.linregress(x,y)

        slope, yint, r = cld_loglinreg[0:3]
        r_squared = r**2.0
        A = np.exp(yint)

        func_form = 'CLD(E, J={J}, Pi={Pi}) = {A}*e^({B}*E), r^2={r_squared}'
        out = func_form.format(J=j, Pi=pi, A=A, B=slope, r_squared=r_squared)

        return out




################################################################################
class ParameterFits:
     """
     Inputs adjusted or unadjusted level density arrays and fits level density
     arrays to rho(E, J, Pi) functions. Determines level density parameters to
     accomplish the function using the various LD models
     """
################################################################################
