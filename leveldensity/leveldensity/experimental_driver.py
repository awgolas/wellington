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
################################################################################

def run():

    target = '52Cr'
    
