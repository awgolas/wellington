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
from scipy import integrate

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


    cld_data = {}
    cld_estimate = {}
    for jpi in Jpi:
        if jpi[1] == 0: continue
        if jpi[0] == -1: continue

        lda = experimental_fit.LevelDensityAnalyzer(inputdict, jpi)
        if len(lda.exp_cld) < 10: continue

        x,y = lda.cld_hist
        mc = Math(x_data=x, y_data=y)
        
        cld_array = mc.smoothing_2d(window=5, mode='exponential')
        cld_est = cld_array[0]
        #param_array = cld_array[1]

        #if np.shape(rho_est) == (2,0):
        #    continue

        #print(np.shape(rho_est), 'SHAPE OF ESTIMATED RHO')
        #print(rho_est)

        #cld_est = integrate.cumtrapz(rho_est[1, :], rho_est[0, :])
        #rho_est_array = np.asarray([rho_est[0, 1:], cld_est])
    #    print(np.shape(cld_est))
        label = jpi
        cld_data[label] = lda.cld_hist
        cld_estimate[label] = cld_est
        #print(np.shape(cld_estimate[label]))
        #y,x = cld_yx
        print(label)
        #plt.step(x,y, label=label)
    #keys = [(2.0, 1.0), (3.0, -1.0), (6.0, -1.0), 'total']
#print(cld_xy)
    fig, ax = plt.subplots(nrows=1, ncols=9,
            figsize=(50,5))
    for n, key in enumerate(cld_estimate.keys()):
        x_data, y_data = cld_data[key]
        x_est = cld_est[0, :]
        y_est = cld_est[1, :]
    
        ax[n].plot(x_est,y_est, label='estimate {}'.format(key))

    
        ax[n].step(x_data, y_data, label='actual {}'.format(key))
        ax[n].set_xlabel('Excitation Energy (MeV)')
        ax[n].set_ylabel('Cumulative Level Distribution')
        ax[n].legend()
    #    if key is not 'total':
    #        jeepee = 'J' + str(key[0]) + 'pi' + str(key[1])
    #    else:
    #        jeepee = key
    fig.savefig('/home/agolas/Pictures/ctfspins.png')
    plt.close()


################################################################################
if __name__ == '__main__':
    run()
