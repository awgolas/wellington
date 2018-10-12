################################################################################
# Author: Alec Golas                                                           #
# Date: October 4th, 2018                                                      #
# Title: experimental_fit.py                                                   #
# Fits experimental level densities to curves                                  #
################################################################################

import os

import numpy as np
import sys
from utilities import Math
from library import Loader

################################################################################
class PostLevelDensityLoad:
    """
    Pulls the levels for a nucleus from the RIPL or OSLO gamma emissions
    libraries and sorts them into CLD(E, J) and CLD(E) data sets
    """

    def __init__(self, parameters):
        self.parameters = Loader(parameters).parameters

    def get_data_parameters(self, source='RIPL'):

        target = self.parameters['target_label']
        parity = self.parameters['parity']
        spin  = self.parameters['spin']



        #data source
        #spin, parity, energy-range

    def ripl_data(self):
        #Open file with z{charge number}.dat formatting
        #Parse for target nucleus symbol
        #Ignore gamma emissions data
        #output data array with [level index, energy, spin, parity]
        #Also get D0 from other different file source

    def oslo_data(self):
        #Do check to make sure nuclei is in the oslo folder
        #Open file with {nucleisymbol}.dat
        #That's basically it, they're preformatted
        #Also calculate CLD if it is requested
        #output data array as [level index, energy, level density]

    def ripl_sort(self, **kwargs):

        #possibly accepted keywords:
            #spin
            #parity
        #calculates cumulative level density
        #calculates rho(E) from the derivative of the cumulative level density
        #outputs array as [excitation energy, spin, parity, rho]

################################################################################
class LevelDensityAnalyzer(PostLevelDensityLoad):

    def __init__(self, levels):
        self.levels = levels
        self.ex_energy = ex_energy

    @property
    def num_levels(self):
        return len(self.levels)

    @property
    def cumulative_levels(self):
        cld = Math().integral_linear(self.ex_energy, self.levels)
        return cld

    @property
    def spacings(self):
        return 1/self.levels

    @property
    def num_spacings(self):
        return len(self.spacings)

    @property
    def ave_spacing(self):
        return numpy.average(self.spacings)

    @property
    def mean_spacing(self):
        return numpy.mean(self.spacings)

    @property
    def stddev_spacing(self):
        return numpy.std(self.spacings)

    @property
    def var_spacing(self):
        return numpy.var(self.spacings)

    @property
    def normalized_spacings(self):
        result = []
        for i in range(self.num_levels-1):
            D = self.levels[i+1]-self.levels[i]
            E = 0.5*(self.levels[i+1]+self.levels[i])
            result.append(D/self.mean_spacing_at_E(E))
        return result

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
