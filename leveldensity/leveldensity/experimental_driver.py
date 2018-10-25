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

sys.path.append('/home/agolas/fudge')
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

from utilities import Math, CurveFitting
import experimental_fit
import level_parameters
################################################################################

def run():

    target = '51Cr'
    inputdict = {'target' : target}

    pldl = experimental_fit.PostLevelDensityLoad(inputdict)
    #Jpi = pldl.Jpi
    Jpi = ['total']
    #print(Jpi)
    #rd = pldl.exp_cld[(1.5, -1.0)]
    #print(rd)


    cld_data = {}
    cld_estimate = {}
    param_estimate = {}

    for jpi in Jpi:
        if jpi[1] == 0: continue
        if jpi[0] == -1: continue

        lda = experimental_fit.LevelDensityAnalyzer(inputdict, jpi)
        if len(lda.exp_cld) < 10: continue

        x,y = lda.cld_hist[:, :150]
        
        pe = experimental_fit.ParameterEstimates(inputdict, jpi)
        A = pe.spin_parity_rho
        B = 1.0/pe.temperature
        C = 1.80
        globalparams = [A, B, C]
        mc = CurveFitting(x_data=x, y_data=y, global_params=globalparams)
        cld_array = mc.smoothing_2d(window=5, mode='exponential')
        cld_est = cld_array[0]
        param_array = cld_array[1]

        A_est = param_array[:,0]
        B_est = param_array[:,1]
        C_est = param_array[:,2]

        spin_dep_calc = A_est + globalparams[0]
        temp_calc = 1/(B_est + globalparams[1])
        cutoff_calc = C_est + globalparams[2]
        #print(cutoff_calc)

        label = jpi
        cld_data[label] = lda.cld_hist
        cld_estimate[label] = cld_est
        param_estimate[label] = [spin_dep_calc, temp_calc, cutoff_calc]

    fig, ax = plt.subplots(nrows=1, ncols=2,
            figsize=(12,5))
    for n, key in enumerate(cld_estimate.keys()):
        x_data, y_data = cld_data[key][:, :150]
        x_est = cld_est[0, :]
        y_est = cld_est[1, :]

        spin_dep, temp, cutoff = param_estimate[key]

        ax[0].step(x_data, y_data, label='actual {}'.format(key))
        ax[0].semilogy(x_est, y_est, label='Calculated w/ Energy-Dependent Parameters')

        ax[1].plot(x_data, spin_dep, label='Spin-Dependence')
        ax[1].plot(x_data, temp, label='Temperature')
        ax[1].plot(x_data, cutoff, label='Cutoff Energy')

        ax[0].set_xlabel('Excitation Energy (MeV)')
        ax[1].set_xlabel('Excitation Energy (MeV)')

        ax[0].set_ylabel('Cumulative Level Distribution')
        ax[1].set_ylabel('Parameter Value (MeV)')

        ax[0].legend()
        ax[1].legend()

        ax[0].grid()
        ax[1].grid()
    #    if key is not 'total':
    #        jeepee = 'J' + str(key[0]) + 'pi' + str(key[1])
    #    else:
    #        jeepee = key
    fig.savefig('/home/agolas/Pictures/ctfspins.png')
    plt.close()

def histogram():

    target = '53Cr'
    inputdict = {'target' : target}

    ripl_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='RIPL')
    hfb_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='HFB')

    err_jpis = ['total', (-1.0, 0.0), (7.0, 0.0), (2.0, 0.0), (9.0, 0.0),
            (10.0, 0.0), (11.0, 0.0)]
    jpis = [i for i in ripl_lda.Jpi if i not in err_jpis]

    hfb_vals = []
    ripl_vals = []

    for jpi in jpis:
        rval = ripl_lda.cld_extract(ex_energy=14.0, j_pi=jpi)
        hval = hfb_lda.cld_extract(ex_energy=14.0, j_pi=jpi)

        ripl_vals.append(rval)
        hfb_vals.append(hval)

    hfb_norm = Math(array=hfb_vals).normalize()
    ripl_norm = Math(array=ripl_vals).normalize()

    jxpis = []

    for js in jpis:
        jxpi = js[0]*js[1]
        jxpis.append(jxpi)


    title_form = """Normalized Distribution of $^{54}Cr$ at U=14.0 MeV """# \n


    plt.bar(jxpis, hfb_norm, align='edge', width=-0.4, label='HFB')
    plt.bar(jxpis, ripl_norm, align='edge', width=0.4, label='RIPL')


    plt.xlabel('$J*\pi$')
    plt.ylabel('$CLD/NCUMUL$')
    plt.title(title_form)
    plt.legend()

    plt.grid()
    plt.savefig('/home/agolas/Pictures/cr54U14.png')
    plt.close()


def knownvsunknown():

    target = '51Cr'
    inputdict = {'target' : target}

    pldl = experimental_fit.PostLevelDensityLoad(inputdict)
    Jpi = pldl.Jpi

    jpitot = 'total'

    lda = experimental_fit.LevelDensityAnalyzer(inputdict, jpitot)
    xtot,ytot = lda.cld_hist

    ex_eng = known_cld = []


    for jpi in Jpi:
        if jpi == 'total':
            continue
        if jpi[0] == -1:
            continue
        if jpi[1] == 0.:
            continue
        lda = experimental_fit.LevelDensityAnalyzer(inputdict, jpi)

        xjpi, yjpi = lda.cld_hist
        for val in xjpi:
            ex_eng.append(val)
        
    num_lvls = len(ex_eng)
    known_cld = range(0, num_lvls+1)
    plt.step(xtot, ytot, label='All levels in RIPL')
    plt.step(sorted(ex_eng), known_cld[1:], label='All levels with defined J & $\pi$')

    plt.grid()
    plt.legend()
    plt.xlabel('Excitation Energy (MeV)')
    plt.ylabel('CLD')
    plt.savefig('/home/agolas/Pictures/known.png')
    


################################################################################
if __name__ == '__main__':
    histogram()
