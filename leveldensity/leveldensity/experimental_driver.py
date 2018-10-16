################################################################################
# Author: Alec Golas                                                           #
# Date: October 15, 2018                                                       #
# Title: experimental_driver.py                                                #
# Drives the experimental fitting module                                       #
################################################################################
from __future__ import print_function

import os
import sys
import math

import numpy as np
import json

from utilities import Math
import experimental_fit
import level_parameters
################################################################################

def run():

    target = '51Cr'
    inputdict = {'target' : target}

    pldl = experimental_fit.PostLevelDensityLoad(inputdict)
    Jpi = pldl.Jpi
    #print(Jpi)
    #rd = pldl.exp_cld[(1.5, -1.0)]
    #print(rd)

    calc_levels = 'literally nothing rn'

    for jpi in Jpi:
        lda = experimental_fit.LevelDensityAnalyzer(inputdict, jpi, calc_levels)

        cep = lda.cld_equation_parameters
        print(cep)

################################################################################
if __name__ == '__main__':
    run()


