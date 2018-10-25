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
import scipy.integrate
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

    def __init__(self, inputdict, source):
        Parameters.__init__(self, inputdict)
        self.exp_cld = self.get_data(source)
        self.source = source
        

    @property
    def RIPL_dir(self):
        return os.path.join('/home/agolas', 'empire', 'RIPL')

    @property
    def levels_dir(self):
        return os.path.join(self.RIPL_dir, 'levels')

    @property
    def HFB_dir(self):
        return os.path.join(self.RIPL_dir, 'densities','total',
                'level-densities-hfb/')


    def get_data(self, source='RIPL'):

        nucleus = (self.compound_label, self.mass_number, self.num_protons)
        parity = self.pi
        spin  = self.spin

        if source == 'RIPL':
            data = self.ripl_data
        elif source == 'Oslo':
            data = self.oslo_data
        elif source == 'HFB':
            data = self.hfb_data
        else:
            raise KeyError("Level Density Data Source must be either 'RIPL',\
                    'HFB', or 'Oslo'")

        return data


    @property
    def ripl_data(self):
            
        leveldir = self.levels_dir
        
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
    def hfb_spin_vals(self):
        j_vals = np.linspace(0,49, num=50)
        if self.mass_number % 2 != 0:
            j_vals = j_vals + 0.5
        return j_vals


    @property
    def hfb_table(self):

        Z_fill = str(self.num_protons).rjust(3)
        A_fill = str(self.mass_number).rjust(3)
        negParityToken = "Z={} A={}: Negative-parity".format(Z_fill, A_fill)
        posParityToken = "Z={} A={}: Positive-Parity".format(Z_fill, A_fill)

        file_name = 'z{:03d}.tab'.format(self.num_protons)
        hfb_path = self.HFB_dir + file_name
        f = open(hfb_path, 'r').read()

        hfb_table = {}
        hfb_table['U'] = []
        hfb_table['T'] = []

        hfb_table[('total', 1)] = []
        hfb_table[('total', -1)] = []
        
        lines = f.split('\n')

        j_vals = self.hfb_spin_vals

        #Initialize the loop so that the lists can just be appended
        for pi in [-1,1]:
            for J in j_vals:
                hfb_table[(J,pi)] = []

        # Read in the data tables
        while lines:
            line = lines.pop(0)

            # Positive parity entries come first
            if posParityToken in line:
                #theTable[1] = []
                lines.pop(0)  # the "****" line
                line = lines.pop(0)  # The column headings

                while line.strip() != '':
                    line = lines.pop(0)
                    sline = line.split()

                    if len(sline) == 0: continue
                    if sline[0] == 'U[MeV]': continue

                    hfb_table['U'].append(float(sline[0]))
                    hfb_table['T'].append(float(sline[1]))
                    hfb_table[('total', 1)].append(float(sline[2]))

                    jcols = sline[5:]

                    for i, j in enumerate(jcols):

                        hfb_table[(j_vals[i],1)].append(float(j))


            # Negative parity entries come second, we're done once we've read these
            if negParityToken in line:

                lines.pop(0)  # the "****" line
                
                while line.strip() != '':
                    line = lines.pop(0)
                    sline = line.split()

                    if len(sline) == 0: continue
                    if sline[0] == 'U[MeV]': continue

                    hfb_table[('total', -1)].append(float(sline[2]))

                    jcols = sline[5:] 
                    for i, j in enumerate(jcols):
                        hfb_table[(j_vals[i],-1)].append(float(j))
                break

        return hfb_table

    @property
    def hfb_temp(self):

        mapped_temp = np.asarray([self.hfb_table['U'], self.hfb_table['T']])
        return mapped_temp

    @property
    def hfb_eff_energy(self):
        return self.hfb_table['U']

    @property
    def hfb_data(self):

        hfb_table = self.hfb_table
        j_vals = self.hfb_spin_vals
        
        hfb_cld = {}
        hfb_rho = {}

        interp_j = np.linspace(0.0,49,num=99)

        for pi in [-1,1]:
            for j in j_vals:

                if (j,pi) in hfb_table.keys():

                    hfb_rho[(j,pi)] = hfb_table[(j,pi)]
                else:

                    rho_l = hfb_table[(j-0.5,pi)]
                    rho_r = hfb_table[(j+0.5,pi)]
                    rho_ave = np.sum([rho_l, rho_r], axis=0)/2

                    hfb_rho[(j,pi)] = rho_ave

        x = hfb_table['U']
        for pi in [-1,1]:
            for j in j_vals:
                y = hfb_rho[(j,pi)]

                cld = scipy.integrate.cumtrapz(y, x, initial=0)
                hfb_cld[(j,pi)] = cld


        pneg =  hfb_table[('total',-1)]
        ppos =  hfb_table[('total', 1)]
        total_cld = np.sum([pneg, ppos], axis=0)
        hfb_cld['total'] = total_cld

        return hfb_cld


    @property
    def Jpi(self):

        cld_hist = self.get_data(source=self.source)

        j_pis = list(cld_hist.keys())
        return j_pis

    def oslo_data(self, nucleus):
        #Do check to make sure nuclei is in the oslo folder
        #Open file with {nucleisymbol}.dat
        #That's basically it, they're preformatted
        #Also calculate CLD if it is requested
        #output data array as [level index, energy, level density]
        return 'Error my guy'

################################################################################
class LevelDensityModelEstimates(PostLevelDensityLoad):
    """
    Inputs adjusted or unadjusted level density arrays and fits level density
    arrays to rho(E, J, Pi) functions. Determines level density parameters to
    accomplish the function using the various LD models. 
    """

    def __init__(self, inputdict, source):
        PostLevelDensityLoad.__init__(self, inputdict, source)
        pldl = PostLevelDensityLoad(inputdict, source)

    @property
    def eff_energy_correction(self):
        return 0.173015

    @property
    def spin_parity_rho(self):

        jpi = self.Jpi
        sigma2 = self.global_spin_cutoff

        spin_parity = {}

        for j_pi in jpi:

            if j_pi =='total':
                rho_jp = 1/math.sqrt(2*sigma2)

            else:
                j = j_pi[0]
                rho_jp = np.exp(-1.0*(j+0.5)**(2.0)/(2.0*sigma2))\
                *(0.5)*(2.0*j+1)/((8.0*sigma2**(1.5))**(0.5))

            spin_parity[j_pi] = rho_jp

        return spin_parity

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

class LevelDensityAnalyzer(LevelDensityModelEstimates):

    def __init__(self, inputdict, source):
        LevelDensityModelEstimates.__init__(self, inputdict, source)
        pldl = PostLevelDensityLoad(inputdict, source)

    @property
    def cld_xy(self):

        cld_mapped_jpis = {}

        for jpi in self.Jpi:

            exp_cld = self.exp_cld[jpi]

            cld_energy = np.asarray(exp_cld)

            if cld_energy[0] == 0:
                cld_energy = cld_energy[1:]

            if self.source == 'RIPL':
                num_levels = len(cld_energy)
                cld = np.arange(1, num_levels+1, step=1)
                mapped = np.array([cld_energy, cld])

            elif self.source =='HFB':
                eff_energy = self.hfb_eff_energy
                cld = self.exp_cld[jpi]
                mapped = np.array([eff_energy, cld])

            cld_mapped_jpis[jpi] = mapped


        return cld_mapped_jpis

    def cld_extract(self, ex_energy, j_pi='total'):

        energy, cld = self.cld_xy[j_pi]
        index = Math(array=energy, value=ex_energy).nearest_value_index()

        extracted_cld = cld[index]

        return extracted_cld


    @property
    def rho_two_equation(self):

        ex_energy, cld  = self.cld_hist

        rho_energy, rho = Math(x_data=ex_energy, y_data=cld).dfdx_1d()

        e_m = self.matching_energy


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

################################################################################
