################################################################################
# Author: Alec Golas                                                           #
# Date: October 4th, 2018                                                      #
# Title: experimental_fit.py                                                   #
# Fits experimental level densities to curves                                  #
################################################################################

import os
import sys
import re
import math

import numpy as np
from scipy.optimize import curve_fit

from utilities import Math
from library import Parameters
from level_parameters import CompositeGilbertCameronParameters as CGCP
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

        jpi_col = ['total',]
        cld_hist = {}
        cld_hist['total'] = []
        for line in formatted_nucleus_data:
            sline = line.split()

            index = sline[0]
            eff_energy = float(sline[1])
            J = float(sline[2])
            pi = float(sline[3])

            Jpi = (J,pi)


            cld_hist['total'].append(eff_energy)
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
class ParameterEstimates(PostLevelDensityLoad):
    """
    Inputs adjusted or unadjusted level density arrays and fits level density
    arrays to rho(E, J, Pi) functions. Determines level density parameters to
    accomplish the function using the various LD models
    """

    def __init__(self, inputdict, JPi):
        PostLevelDensityLoad.__init__(self, inputdict)
        pldl = PostLevelDensityLoad(inputdict)
        self.JPi = JPi

    @property
    def eff_energy_correction(self):
        return 0.173015

    @property
    def spin_parity_rho(self):

        sigma2 = self.global_spin_cutoff

        if self.JPi =='total':
            rho_jp = 1/math.sqrt(2*sigma2)

        else:
            j = self.JPi[0]
            rho_jp = np.exp(-1.0*(j+0.5)**(2.0)/(2.0*sigma2))\
            *(0.5)*(2.0*j+1)/((8.0*sigma2**(1.5))**(0.5))

        return rho_jp

    @property
    def atilda(self):

        A = self.mass_number
        atilda = 0.154*A + 6.3e-5*A**2
        return atilda

    @property
    def gamma(self):
        gam0 = 0.41029
        A = self.mass_number

        g = gam0/A**0.333
        return g

    @property
    def temperature(self):
        A = self.mass_number
        g = self.gamma
        delta_W = self.shell_correction

        temp = -0.22 + 9.4/math.sqrt(A*(1.0 + g*delta_W))
        return temp


#################################################################################

class LevelDensityAnalyzer(ParameterEstimates):

    def __init__(self, inputdict, JPi):
        ParameterEstimates.__init__(self, inputdict, JPi)
        pldl = PostLevelDensityLoad(inputdict)
        self.exp_cld = pldl.exp_cld[JPi]

    @property
    def cld_hist(self):

        cld_energy = np.asarray(self.exp_cld)

        if cld_energy[0] == 0:
            cld_energy = cld_energy[1:]

        num_levels = len(cld_energy)
        index = np.arange(1, num_levels+1, step=1)

        histogram = np.array([cld_energy, index])
        return histogram

    @property
    def rho_two_equation(self):

        ex_energy, cld  = self.cld_hist

        rho_energy, rho = Math().dfdx_1d(ex_energy, cld)

        e_m = self.matching_energy

        fgf = self.fermi_gas_rho_form
        ctf = self.constant_temp_rho_form

        if np.min(rho_energy) > e_m:
            rho_estimate = self.rho_estimation(fgf, rho_energy, rho)
        elif np.max(rho_energy) < e_m:
            rho_estimate = self.rho_estimation(ctf, rho_energy, rho)
        else:
            cutoff = Math().nearest_value_index(rho_energy, e_m)

            low_energy = rho_energy[:cutoff]
            high_energy = rho_energy[cutoff:]
            low_rho = rho[:cutoff]
            high_rho = rho[cutoff:]
            try:
                rho_low_est = self.rho_estimation(ctf, low_energy, low_rho)
            except:
                rho_low_est = np.asarray([[], []])
            try:
                rho_high_est = self.rho_estimation(fgf, high_energy, high_rho)
            except:
                rho_high_est = np.asarray([[], []])

            rho_estimate = np.hstack((rho_low_est, rho_high_est))

        return rho_estimate


    def cld_smoother(self):
        ex_energy , index = self.cld_hist

        mc = Math(x_data=ex_energy, y_data=index)
        cld_curve = mc.smoothing_2d(window=9, mode='exponential')

        return cld_curve


    def rho_estimation(self, curveform, x_data, y_data):

        jpi = self.JPi

        rho_jpi = self.spin_parity_rho


        if curveform == self.constant_temp_rho_form:

            temp = self.temperature
            p0 = [rho_jpi, 1, temp, 1]

        elif curveform == self.fermi_gas_rho_form:

            atilda = self.temperature
            gamma = self.gamma
            correction = self.eff_energy_correction
            p0 = [rho_jpi, atilda, gamma, correction, 1]
            
        popt, pcov = curve_fit(curveform, x_data, y_data, p0=p0,maxfev=10000000)

        y_estimate = curveform(x_data, *popt)
        rho_estimate = np.asarray([x_data, y_estimate])

        return rho_estimate

    def constant_temp_rho_form(self, E, spindep, cutoff, temp, const):

        return spindep/temp*np.exp((E-cutoff)/temp) + const

    def fermi_gas_rho_form(self, E, spindep, atilda, gamma, correction,  const):
        # A --> spin dependence
        # B --> atilda
        # C --> gamma
        # D --> pairing energy correction
        # F --> fitting constant

        delta_W = self.shell_correction
        delta = self.delta

        U = E - delta - correction

        a_param = atilda * (1 + (1 - np.exp(gamma*U)*delta_W/U))

        func = spindep/(12.0*a_param**0.25*U**1.25)\
                *np.exp(2.0*np.sqrt(a_param*U)) + const
        return func

################################################################################
