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


################################################################################
class LevelDensityLoad:
    """ Pulls the levels for a nucleus from the RIPL or OSLO gamma emissions
    libraries and sorts them into CLD(E, J) and CLD(E) data sets"""
    def get_data_parameters(self):

        #target nuclei
        #data source
        #spin, parity, energy-range

    def ripl_data(self):
        #Open file with z{charge number}.dat formatting
        #Parse for target nucleus symbol
        #Ignore gamma emissions data
        #output data array with [level index, energy, spin, parity]

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
class LevelDensityAdjustments:
    """Inputs formatted array data that has been loaded from LevelDensityLoad and
    adjusts the level densities to account for missing levels. This only accepts
    formatted RIPL data, as Oslo-method data should theoretically be completely
    defined.

    The getFractionMissing2SpacingCorrelation and get2get2SpacingCorrelation
    functions were modified from the fudge-4.2.3 package, BNL.restools"""



    def getFractionMissing2SpacingCorrelation(self):
        """
        See G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using Neutron Spacings",
        IAEA NDS Report INDC(NDS)-0561 (2009)
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
        """
        Get the spacing-spacing correlation function.

            ..math::
                D_i=E_{i+1}-E_i

            ..math::

                \rho(D_i,D_{i+1}) = \frac{\sum_i(D_i-\overline{D})(D_{i+1}-\overline{D})}{\left[\sum_i(D_i-\overline{D})^2\sum_j(D_{j+1}-\overline{D})^2\right]^{1/2}}

        If you stare at this for a while, it is just the level-level level spacing correlation coefficient.

        For GOE level spacings, should see :math:`\rho=-0.27`

        :returns : value and variance in a tuple
        """
        diff_spacings=[]
        for i in range(self.num_levels-1):
            aveE=0.5*(self.levels[i+1]+self.levels[i])
            diff_spacings.append(self.spacings[i]-self.mean_spacing_at_E(aveE))
        if len([x for x in diff_spacings if x>epsD])<2:
            raise ValueError("Level-level spacing correlation undefined, is this a picket fence?")
        correlation_matrix = numpy.corrcoef(numpy.array([[diff_spacings[i], diff_spacings[i+1]] for i in range(self.num_spacings-1)]).T)
        return correlation_matrix[0,1]

    pass
################################################################################
class ParameterFits:
    """Inputs adjusted or unadjusted level density arrays and fits level density
    arrays to rho(E, J, Pi) functions. Determines level density parameters to
    accomplish the function using the various LD models"""

    pass
################################################################################
