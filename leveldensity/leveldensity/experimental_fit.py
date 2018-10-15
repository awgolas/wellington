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
        

    @property
    def RIPL_dir(self):
        return os.path.join('~', 'empire', 'RIPL', 'levels')

    def get_data(self, source='RIPL'):

        nucleus = (self.compound_label, self.mass_number, self.num_protons)
        parity = self.pi
        spin  = self.spin

        if source == 'RIPL':
            data = self.ripl_data(nucleus, parity, spin)
        elif source == 'Oslo':
            data = self.oslo_data(nucleus)


        #data source
        #spin, parity, energy-range

    def ripl_data(self, nucleus, parity, spin):
            
        leveldir = self.RIPL_dir
        
        label = nucleus[0]
        r = re.compile("([0-9]+)([a-zA-Z]+)")
        m = r.match(label)
        end_z = str(int(m.group(1)) + 1)
        end_label = end_z + m.group(2)
        

        mass_number = nucleus[1]
        charge_number = nucleus[2]

        file_name = 'z{:03d}.dat'.format(charge_number)
        with open(leveldir + file_name, 'r') as f:
            data_file = list(f)

        target = False
        for n, line in enumerate(data_file):
            
            if label in line:
                target = True
                init = n
            
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
        cld_hist = {}
        for line in formatted_nucleus_data:
            sline = line.split()
            
            eff_energy = sline[1]
            J = sline[2]
            pi = sline[3]

            Jpi = (J,pi)
            if Jpi not in jpi_col:
                jpi_col.append(Jpi)
                cld_hist[Jpi] = []

            cld_hist[Jpi].append(eff_energy)
         
        return cld_hist



        #Open file with z{charge number}.dat formatting
        #Parse for target nucleus symbol
        #Ignore gamma emissions data
        #output data array with [level index, energy, spin, parity]
        #Also get D0 from other different file source

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

    def __init__(self, inputval, calc_levels):
        PostLevelDensityLoad.__init__(self, inputval)
        pldl = PostLevelDensityLoad(inputval)
        self.exp_levels = pldl.exp_levels
        self.calc_levels = calc_levels

    @property
    def num_exp_levels(self):
        return len(self.exp_levels)

    @property
    def cumulative_levels(self):
        cld = Math().integral_linear(self.excitation_energy, self.calc_levels)
        return cld

    @property
    def exp_spacings(self):
        return 1/self.exp_levels

    @property
    def num_exp_spacings(self):
        return len(self.spacings)

    @property
    def mean_spacing(self):
        return numpy.mean(self.exp_spacings)

    @property
    def stddev_exp_spacing(self):
        return numpy.std(self.exp_spacings)

    @property
    def var_exp_spacing(self):
        return numpy.var(self.exp_spacings)




################################################################################
class LevelDensityAdjustments(LevelDensityAnalyzer):
    """
    Inputs formatted array data that has been loaded from LevelDensityLoad and
    adjusts the level densities to account for missing levels. This only accepts
    formatted RIPL data, as Oslo-method data should theoretically be completely
    defined.

    The getFractionMissing2SpacingCorrelation and get2get2SpacingCorrelation
    functions were modified from the fudge-4.2.3 package, BNL.restools
    """

    def getFractionMissing2SpacingCorrelation(self):
        """
        See G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using
        Neutron Spacings", IAEA NDS Report INDC(NDS)-0561 (2009)
        """
        rho=self.get2SpacingCorrelation()
        emean=numpy.array([-0.251, 0.428])
        ecov=numpy.array([[2.67e-4, -7.44e-4], [-7.44e-4, 3.22e-3]])
        dfde0=-1.0/emean[1]
        dfde1=-(rho-emean[0])/emean[1]/emean[1]
        f=(rho-emean[0])/emean[1]
        df=math.sqrt(dfde0*dfde0*ecov[0][0] + dfde1*dfde1*ecov[1][1] + 2.0*dfde0*dfde1*ecov[0][1])
        return max(min(1.0,f),0.0), df

    def get2SpacingCorrelation(self, epsD=1e-9):

        diff_spacings=[]
        for i in range(self.num_levels-1):

            aveE=0.5*(self.levels[i+1]+self.levels[i])

            diff_spacings.append(self.spacings[i]-self.mean_spacing_at_E(aveE))

        if len([x for x in diff_spacings if x>epsD])<2:

            raise ValueError("Level-level spacing correlation undefined, is this a picket fence?")
        
        correlation_matrix = numpy.corrcoef(numpy.array([[diff_spacings[i], diff_spacings[i+1]] for i in range(self.num_spacings-1)]).T)
   
        return correlation_matrix[0,1]

# ################################################################################
# class ParameterFits:
#     """
#     Inputs adjusted or unadjusted level density arrays and fits level density
#     arrays to rho(E, J, Pi) functions. Determines level density parameters to
#     accomplish the function using the various LD models
#     """
#
#     pass
# ################################################################################
