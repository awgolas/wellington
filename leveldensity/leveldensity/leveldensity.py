################################################################################
# Author: Alec Golas                                                           #
# Date: September 18th, 2018                                                   #
# Title: leveldensity.py                                                       #
# Determines cr52 Level Densities for plotting using different models          #
################################################################################
from __future__ import print_function
import matplotlib
matplotlib.use('agg')

import os
import sys
import numpy as np
import math
import json
import matplotlib.pyplot as plt
from utilities import Math
from library import Loader

import level_parameters

def run():

    excitation_energy = np.linspace(3.85,10,num=50)
    parity = 1
    target = '52Cr'
    projectile = 'neutron'
    spin = [0.0, 1.0, 1.5, 2.0]
    fig, ax = plt.subplots(nrows=1, ncols=1)
    for j in spin:
        input = {'target' : target,
                 'spin'   : j,
                 'parity' : parity,
                 'excitation_energy' : excitation_energy}


        fgm = level_parameters.BackShiftedFermiGasParameters(input)
        fg_energy = fgm.eff_energy
        fgm_d = 1/fgm.leveldensity

        gcm = level_parameters.CompositeGilbertCameronParameters(input)
        gc_energy = gcm.eff_energy
        gcm_d = 1/gcm.leveldensity

        ax.semilogy(fg_energy, fgm_d, label='FG 1/rho(u, $J$={}, $\pi_0$)'.format(j))
        ax.semilogy(gc_energy, gcm_d, label='GC 1/rho(u, $J$={}, $\pi_0$)'.format(j))

    save_dir = '/home/agolas/Documents/levden'

    d0_array = 0.032*np.ones(np.shape(fgm_d))
    d0_error = 0.0035*np.ones(np.shape(fgm_d))
    eff_sep_energy = gcm.separation_energy - gcm.delta
    bn = [eff_sep_energy]
    bn_err =
    d0 = 0.032

    ax.errorbar(fg_energy, d0_array, yerr=d0_error, label='RIPL $D_0$')
    ax.semilogy(bn_array, [0.0,0.35], label='Neutron Separation Energy')
    ax.legend(loc=0)
    ax.grid()
    ax.set_xlabel('Effective Excitation Energy $MeV$')
    ax.set_ylabel('Level Spacing $MeV$')
    ax.set_title('Level Spacing of $^{52}Cr$')
    fig.savefig(save_dir + 'mls.png')

if __name__ == "__main__":
    run()
################################################################################
######################## End of leveldensity.py ################################
