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

    zvv_file = """4.19   2.21655e+01
4.53   3.59912e+01
4.88   5.76440e+01
5.23   7.65317e+01
5.58   9.38064e+01
5.93   1.13063e+02
6.28   1.35094e+02
6.63   1.61073e+02
6.98   1.92428e+02
7.33   2.30762e+02
7.67   2.77853e+02
8.02   3.35690e+02
8.37   4.06536e+02
8.72   4.92999e+02"""

    lines = zvv_file.split('\n')
    ldy = []
    ldx = []
    for line in lines:
        sline = line.split('   ')
        ld = float(sline[1])
        ex = float(sline[0])
        ldx.append(ex)
        ldy.append(ld)



    excitation_energy = np.linspace(1.85,10,num=50)
    parity = 1
    target = '52Cr'
    spin = [1.0]#0.5, 1.0, 1.5, 2.0]
    for j in spin:
        input = {'target' : target,
                 'spin'   : j,
                 'parity' : parity,
                 'excitation_energy' : excitation_energy}

        #p = plt.scatter(ldx, ldy, label='EMPIRE EGSM rho(u)')
        #p = plt.scatter(fg_energy, fgm, label='FG rho(u,$I_0$, $\pi_0$)')
        #p = plt.scatter(gc_energy, gcm, label='GC rho(u,$I_0$, $\pi_0$)')
        #p = plt.scatter(gc_energy, gcm_rhoe, label='GC rho(u)')
    #plt.legend(loc=0)
    #plt.grid()
    #plt.xlabel('Effective Excitation Energy $MeV$')
    #plt.ylabel('Level Density $MeV^{-1}$')
    #p.figure.savefig('rhojpi.png')

if __name__ == "__main__":
    run()
################################################################################
######################## End of leveldensity.py ################################
