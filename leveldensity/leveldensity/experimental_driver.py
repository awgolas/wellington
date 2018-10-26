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

def histogram(target):

    inputdict = {'target' : target}

    ripl_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='RIPL')
    hfb_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='HFB')

    A = ripl_lda.mass_number
    Z = ripl_lda.num_protons
    N = A - Z

    jpis = []
    for Js in ripl_lda.Jpi:
        if Js[0] == -1: continue
        if Js[1] == 0: continue
        if Js == 'total': continue
        jpis.append(Js)

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

    error_per_jpi = np.abs(hfb_norm - ripl_norm)
    maximum_error = np.max(error_per_jpi)
    cumulative_error = np.sum(error_per_jpi)

    return (Z, N, A, cumulative_error, maximum_error)

#    for js in jpis:
#        jxpi = js[0]*js[1]
#        jxpis.append(jxpi)
#
#
#    title_form = """Normalized Distribution of $^{54}Cr$ at U=14.0 MeV """# \n
#
#
#    plt.bar(jxpis, hfb_norm, align='edge', width=-0.4, label='HFB')
#    plt.bar(jxpis, ripl_norm, align='edge', width=0.4, label='RIPL')
#
#
#    plt.xlabel('$J*\pi$')
#    plt.ylabel('$CLD/NCUMUL$')
#    plt.title(title_form)
#    plt.legend()
#
#    plt.grid()
#    plt.savefig('/home/agolas/Pictures/cr54U14.png')
#    plt.close()


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

    targets = ['1H', '2H', '3He', '4He', '6Li', '7Li', '9Be', '10B', '11B',
    '12C', '13C', '14N', '15N', '16O', '17O', '18O', '19F', '20Ne', '22Ne',
    '23Na', '24Mg', '25Mg', '26Mg', '27Al', '28Si',  '29Si', '30Si', '31P',
    '32S', '33S', '34S', '35Cl', '36Ar', '36S', '37Cl', '38Ar', '39K', '40Ca',
    '40K', '40Ar', '41K', '42Ca', '43Ca', '44Ca', '45Sc', '46Ti', '46Ca',
    '47Ti', '48Ca', '48Ti', '49Ti', '50Ti', '50V', '50Cr', '51V', '52Cr',
    '53Cr', '54Cr', '54Fe', '55Mn', '56Fe', '57Fe', '58Ni', '58Fe', '59Co',
    '60Ni', '61Ni', '62Ni', '63Cu', '64Ni', '64Zn', '65Cu', '66Zn', '67Zn',
    '68Zn', '69Ga', '70Zn', '70Ge', '71Ga', '72Ge', '73Ge', '74Ge', '75As']

    Z_label = {1: 'H', 
               2: 'He', 
               3: 'Li', 
               4: 'Be', 
               5: 'B', 
               6: 'C',
               7: 'N',
               8: 'O',
               9: 'F',
               10: 'Ne',
               11: 'Na',
               12: 'Mg',
               13: 'Al',
               14: 'Si',
               15: 'P',
               16: 'S',
               17: 'Cl',
               18: 'Ar',
               19: 'K',
               20: 'Ca',
               21: 'Sc',
               22: 'Ti',
               23: 'V',
               24: 'Cr',
               25: 'Mn',
               26: 'Fe',
               27: 'Co',
               28: 'Ni',
               29: 'Cu',
               30: 'Zn',
               31: 'Ga',
               32: 'Ge',
               33: 'As',
               34: 'Se',
               35: 'Br',
               36: 'Kr',
               37: 'Rb',
               38: 'Sr',
               39: 'Y',
               40: 'Zr',
               41: 'Nb',
               42: 'Mo',
               43: 'Tc',
               44: 'Ru',
               45: 'Rh',
               46: 'Pd',
               47: 'Ag',
               48: 'Cd',
               49: 'In',
               50: 'Sn',
               51: 'Sb',
               52: 'Te',
               53: 'I',
               54: 'Xe'}

    N_label = range(0,70)



    cum_error_map = max_error_map = np.zeros((58, 74))

    cum_error_map[:] = np.nan
    max_error_map[:] = np.nan

    cum_error_1d = max_error_1d = np.zeros(78)

    A_x = []

    for prot, name in Z_label.items():
        for N_val in N_label:
            a_label = str(prot+N_val)
            target = a_label + name

            try:
                Z, N, A, cum_error, max_error = histogram(target)
                print('Target: {}  Cum Error: {}  Max Error: {}'.format(target, cum_error, max_error))
            except:
                print('Target {} Failed'.format(target))
                continue

            cum_error_map[Z,N] = cum_error
            max_error_map[Z,N] = max_error

    plt.pcolormesh(cum_error_map)
    plt.grid()
    plt.xlabel('Z')
    plt.ylabel('N')
    plt.title('Cumulative CLD Difference of HFB and RIPL data')
    plt.savefig('/home/agolas/Pictures/cum_error_map.png')
    plt.close()

    plt.pcolormesh(max_error_map)
    plt.grid()
    plt.xlabel('Z')
    plt.ylabel('N')
    plt.title('Maximum CLD Difference of HFB and RIPL data')
    plt.savefig('/home/agolas/Pictures/max_error_map.png')
    plt.close()

#    plt.scatter(A_x, max_error_1d)
#    plt.xlabel('A')
#    plt.ylabel('$|CLD(J, \pi)_{HFB} - CLD(J, \pi)_{RIPL}|_{max}$')
#    plt.title('Maximum CLD Error as a Function of A')
#    plt.savefig('/home/agolas/Pictures/max_error_1D.png')
#    plt.close()

#    plt.scatter(A_x, cum_error_1d)
#    plt.xlabel('A')
#    plt.ylabel('$\Sigma_{M,N}|CLD(J_M,\pi_N)_{HFB} - CLD(J_M, \pi_N)_{RIPL}|$')
#    plt.title('Comulative Error over Js, Pis of each Nucleus as a Function of A')
#    plt.savefig('/home/agolas/Pictures/cum_error_1d.png')




