################################################################################
# Author: Alec Golas                                                           #
# Date: October 4th, 2018                                                      #
# Title: experimental_fit.py                                                   #
# Fits experimental level densities to curves                                  #
################################################################################

import os
import sys
import numpy as np
from utilities import Math

################################################################################
class LevelDensityLoad:
    """ Pulls the levels for a nucleus from the RIPL library and sorts it
    into CLD(E, J) and CLD(E) functions"""
    def get_data(self):

        #target nuclei#
        #spin#
        #energy range#
