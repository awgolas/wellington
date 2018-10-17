################################################################################
# Author: Alec Golas                                                           #
# Date: October 15, 2018                                                       #
# Title: experimental_driver.py                                                #
# Drives the experimental fitting module                                       #
################################################################################
from __future__ import print_function
import matplotlib
matplotlib.use('agg')

import os
import sys
import math

import numpy as np
import matplotlib.pyplot as plt

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


    cep_coeffs = {}
    cld_data = {}
    cld_estimate = {}
    for jpi in Jpi:
        if jpi[1] == 0: continue
        if jpi[0] == -1: continue

        lda = experimental_fit.LevelDensityAnalyzer(inputdict, jpi)
        if len(lda.exp_cld) < 10: continue

        cld_yx = lda.cld_hist
        
        try:
            cld_est = lda.cld_two_equation
        except:
            continue

        #print(cld_est)


        label = jpi
        cld_data[label] = cld_yx
        cld_estimate[label] = cld_est

        #y,x = cld_yx
        #plt.step(x,y, label=label)
    keys = [(2.0, 1.0), (3.0, -1.0), (6.0, -1.0)]
    keys = ['total']
#print(cld_xy)
    fig, ax = plt.subplots(nrows=1, ncols=1,
            figsize=(5,5))
    for n, key in enumerate(keys):
        x_data, y_data = cld_data[key]
        x_est, y_est = cld_estimate[key]
    
        ax.plot(x_est,y_est, label='estimate {}'.format(key))

    
        ax.step(x_data, y_data, label='actual {}'.format(key))
        ax.set_xlabel('Excitation Energy (MeV)')
        ax.set_ylabel('Cumulative Level Distribution')
        ax.legend()
    #    if key is not 'total':
    #        jeepee = 'J' + str(key[0]) + 'pi' + str(key[1])
    #    else:
    #        jeepee = key
    fig.savefig('/home/agolas/Pictures/total.png')


################################################################################
if __name__ == '__main__':
    run()
