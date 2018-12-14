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
import multiprocessing as mp
import json

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

    ripl_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='RIPL')
    Jpi = ripl_lda.Jpi
    #Jpi = ['total']
    #print(Jpi)
    #rd = pldl.exp_cld[(1.5, -1.0)]
    #print(ripl_lda.temperature, ripl_lda.cutoff_energy)


    # cld_data = {}
    # cld_estimate = {}
    # param_estimate = {}
    #print(ripl_lda.delta)

    exp_data = ripl_lda.cld_xy
    n_tot = 0
    temp_array = []
    e0_array = []
    weight_array = []
    Bn = 7.939
    D_0 = 35./1000.
    D_0_inv = 1./D_0
    D_0_err = 0.003
    for jpi in Jpi:
        #if jpi[1] == 0: continue
        if jpi[0] == -1: continue
        if jpi == 'total': continue

        x,y = exp_data[jpi]
        if len(x) < 10: continue
        #x = np.asarray(x) + ripl_lda.delta
        #y = np.hstack(([0.],y))

        #pe = ripl_lda.ParameterEstimates(inputdict)
        A = 1 #ripl_lda.spin_parity_rho[jpi]
        B = 1.0/ripl_lda.temperature
        C = ripl_lda.cutoff_energy
        globalparams = np.asarray([A, B, C])
        mc = CurveFitting(x_data=x, y_data=y, global_params=globalparams)

        A_adj,B_adj,C_adj, D_adj = mc.parameter_determination()
        A_calc = A_adj + globalparams[0]
        B_calc = B_adj + globalparams[1]
        C_calc = C_adj + globalparams[2]
        D = D_adj*100.

        temp_calc = 1./B_calc
        e0_calc = C_calc
        #if jpi[0] == 0:
        #    if jpi[1] == 1:
        n_jpi = len(x)
        n_tot = n_tot + n_jpi

        nom_temp = temp_calc*n_jpi
        nom_e0 = e0_calc*n_jpi
        temp_array.append(nom_temp)
        e0_array.append(nom_e0)

        ex_energy = np.linspace(0.5, Bn+1)
        jpi_cld = np.exp((1./temp_calc)*(x - e0_calc)) + D

        jpi_bn = (1./temp_calc)*np.exp((1./temp_calc)*(ex_energy - e0_calc))
        glob_bn = B*np.exp(B*(ex_energy - C))

        plt.plot(ex_energy, jpi_bn, label='{} Dependent Level Density'.format(str(jpi)))
        #plt.plot(ex_energy, glob_bn, label='RIPL Parameterization')
        #plt.step(x, y, label='Empirical Data')
        plt.errorbar(Bn, D_0_inv, xerr=D_0_err, yerr=D_0_err, fmt='*', label='$D_0^{-1}$')
        plt.grid()
        plt.legend()
        plt.yscale('log')
        plt.xlabel('Excitation Energy')
        plt.ylabel('Level Density (1/MeV)')
        plt.title('Comparison of $D_{0}^{-1}$ to Fitted Parameters with $^{52}Cr$')
        plt.savefig('/home/agolas/Desktop/{}D0.png'.format(str(jpi)))
        plt.close()

        # plt.plot(x, jpi_cld+2, label='{} Dependent CLD'.format(str(jpi)))
        # plt.step(x, y+2, label='Empirical Data')
        # #plt.errorbar(Bn, D_0_inv, xerr=D_0_err, yerr=D_0_err, fmt='*', label='$D_0^{-1}$')
        # plt.grid()
        # plt.legend()
        # plt.yscale('log')
        # plt.xlabel('Excitation Energy')
        # plt.ylabel('Level Density (1/MeV)')
        # plt.title('Comparison of Level Index to Fitted Parameters with $^{52}Cr$')
        # plt.savefig('/home/agolas/Desktop/{}cld.png'.format(str(jpi)))
        # plt.close()

        print("""(J,Pi): {}
                Nlevels: {}
                Temp: {}
                E0: {}
                Spacing at Bn: {} """.format(jpi, n_jpi, temp_calc, e0_calc, jpi_bn))

        #print('GCROT {} {} {}'.format(ripl_lda.num_protons, ripl_lda.mass_number, 1./B_calc))
        #print('GCROE0 {} {} {}'.format(ripl_lda.num_protons, ripl_lda.mass_number, C_calc))

    print("Total Number of Levels:", n_tot)
    #print("Nominally Weighted E0:", nom_e0)
    #print("Nominally Weighted Temperature:", nom_temp)
    e0_array = np.asarray(e0_array)/float(n_tot)
    temp_array = np.asarray(temp_array)/float(n_tot)

    e0_tot = C#np.sum(e0_array)
    temp_tot = np.sum(temp_array)
    print("Weighted E0: {} \nWeighted Temperature {}".format(e0_tot,temp_tot))

    #D_0 = 32./1000.
    #D_0_inv = 1./D_0
    #D_0_err = 3.5/1000.
    #Bn = 7.939
    # # #print(D_0)
    # #
    #ex_energy = np.linspace(1., Bn+2)
    ex_energy = x

    #rho_Bn_global = (1/1.06304)*np.exp((1/1.06304)*(ex_energy-1.30035))
    rho_Bn_adj = (1./temp_tot)*np.exp((1./temp_tot)*(Bn - e0_tot))
    d0_calc = 1./rho_Bn_adj
    print("""Level Spacing at Neutron Separation Energy:
                Actual: {}
                Calculated {}""".format(D_0, d0_calc))
    # #
    # cld_global = np.exp((B)*(ex_energy-C))
    cld_adj = np.exp((ex_energy -e0_tot)/temp_tot)
    #
    #plt.scatter(Bn, D_0_inv, label='Neutron Separation Energy')
    # plt.plot(ex_energy, cld_adj, label='Adjusted Parameterization')
    # plt.step(x, y, label='Empirical Data')
    # plt.grid()
    # plt.legend()
    # plt.yscale('log')
    # plt.xlabel('Excitation Energy')
    # plt.ylabel('Cumulative Level Distribition')
    # plt.title('Comparison of RIPL Parameters to Fitted Parameters with $^{52}$Cr')
    # plt.savefig('/home/agolas/Desktop/test.png')
    # plt.close()
    # #
    # #
    # #
    # # plt.errorbar(Bn, D_0_inv, xerr=D_0_err, fmt='*', color='b', label='$D_0^{-1}$')
    # # plt.plot(ex_energy, rho_Bn_global, label='RIPL parameterization')
    # # plt.plot(ex_energy, rho_Bn_adj, label='Adjusted Parameterization')
    # # plt.grid()
    # # plt.legend()
    # # plt.yscale('log')
    # # plt.title('Level Density at Neutron Separation Energy of $^{52}$Cr')
    # # plt.xlabel('Excitation Energy (MeV)')
    # # plt.ylabel('Level Density (1/MeV)')
    # # plt.savefig('/home/agolas/Desktop/D_0_52.png')
    # # plt.close()

def shmun():

    target = '51Cr'
    inputdict = {'target' : target}

    ripl_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='RIPL')
    #Jpi = ripl_lda.Jpi
    Jpi = ['total']
    #print(Jpi)
    #rd = pldl.exp_cld[(1.5, -1.0)]
    #print(rd)


    # cld_data = {}
    # cld_estimate = {}
    # param_estimate = {}
    #print(ripl_lda.delta)

    exp_data = ripl_lda.cld_xy
    for jpi in Jpi:
        if jpi[1] == 0: continue
        if jpi[0] == -1: continue

        x,y = exp_data[jpi]
        if len(x) < 10: continue

        #pe = ripl_lda.ParameterEstimates(inputdict)
        A = 1 #ripl_lda.spin_parity_rho[jpi]
        B = 1.0/ripl_lda.temperature
        C = ripl_lda.cutoff_energy
        globalparams = np.asarray([A, B, C])

        mc_dub_1 = CurveFitting(x_data=x[:125], y_data=y[:125], global_params=globalparams)
        mc_dub_2 = CurveFitting(x_data=x[125:], y_data=y[125:], global_params=globalparams)
        mc_1 = CurveFitting(x_data=x[:25], y_data=y[:25], global_params=globalparams)
        mc_2 = CurveFitting(x_data=x[25:50], y_data=y[25:50], global_params=globalparams)
        mc_3 = CurveFitting(x_data=x[50:75], y_data=y[50:75], global_params=globalparams)
        mc_4 = CurveFitting(x_data=x[75:100], y_data=y[75:100], global_params=globalparams)
        mc_5 = CurveFitting(x_data=x[100:125], y_data=y[100:125], global_params=globalparams)
        mc_6 = CurveFitting(x_data=x[125:], y_data=y[125:], global_params=globalparams)
        ex_1 = x[:25]
        ex_2 = x[25:50]
        ex_3 = x[50:75]
        ex_4 = x[75:100]
        ex_5 = x[100:125]
        ex_6 = x[125:]
        ex_energy = x
    #
    A_dub_1,B_dub_1,C_dub_1 = mc_dub_1.parameter_determination() + globalparams
    A_dub_2,B_dub_2,C_dub_2 = mc_dub_2.parameter_determination() + globalparams

    A_1_adj, B_1_adj, C_1_adj = mc_1.parameter_determination()
    A_2_adj, B_2_adj, C_2_adj = mc_2.parameter_determination()
    A_3_adj, B_3_adj, C_3_adj = mc_3.parameter_determination()
    A_4_adj, B_4_adj, C_4_adj = mc_4.parameter_determination()
    A_5_adj, B_5_adj, C_5_adj = mc_5.parameter_determination()
    A_6_adj, B_6_adj, C_6_adj = mc_6.parameter_determination()

    A_1 = A_1_adj + globalparams[0]
    B_1 = B_1_adj/50. + globalparams[1]
    C_1 = C_1_adj + globalparams[2]

    A_2 = A_2_adj + globalparams[0]
    B_2 = B_2_adj/50. + globalparams[1]
    C_2 = C_2_adj + globalparams[2]

    A_3 = A_3_adj + globalparams[0]
    B_3 = B_3_adj/50. + globalparams[1]
    C_3 = C_3_adj + globalparams[2]

    A_4 = A_4_adj + globalparams[0]
    B_4 = B_4_adj/50. + globalparams[1]
    C_4 = C_4_adj + globalparams[2]

    A_5 = A_5_adj + globalparams[0]
    B_5 = B_5_adj/50. + globalparams[1]
    C_5 = C_5_adj + globalparams[2]

    A_6 = A_6_adj + globalparams[0]
    B_6 = B_6_adj/50. + globalparams[1]
    C_6 = C_6_adj + globalparams[2]

    # A_2,B_2,C_2 = A_2_adj, B_2_adj, C_2_adj/100. + globalparams
    # A_3,B_3,C_3 = A_3_adj, B_3_adj, C_3_adj/100. + globalparams
    # A_4,B_4,C_4 = A_4_adj, B_4_adj, C_4_adj/100. + globalparams
    # A_5,B_5,C_5 = A_5_adj, B_5_adj, C_5_adj/100. + globalparams
    # A_6,B_6,C_6 = A_6_adj, B_6_adj, C_6_adj/100. + globalparams

    # A_2,B_2,C_2 = mc_2.parameter_determination() + globalparams
    # A_3,B_3,C_3 = mc_3.parameter_determination() + globalparams
    # A_4,B_4,C_4 = mc_4.parameter_determination() + globalparams
    # A_5,B_5,C_5 = mc_5.parameter_determination() + globalparams
    # A_6,B_6,C_6 = mc_6.parameter_determination() + globalparams

    print(B_1, B_2, B_3, B_4, B_5)
    print(C_1, C_2, C_3, C_4, C_5)

    cld_1 = np.exp(B_1*(ex_1 - C_1))
    cld_2 = np.exp(B_2*(ex_2 - C_2))
    cld_3 = np.exp(B_3*(ex_3 - C_3))
    cld_4 = np.exp(B_4*(ex_4 - C_4))
    cld_5 = np.exp(B_5*(ex_5 - C_5))
    cld_6 = np.exp(B_6*(ex_6 - C_6))
    cld_tot = np.hstack((np.hstack((np.hstack((cld_1,cld_2)),np.hstack((cld_3,cld_4)))), np.hstack((cld_5,cld_6))))

    cld_global = np.exp(B*(ex_energy-C-2.9))
    #plt.step(x[:125], y[:125], label='Empirical Data')
    plt.plot(ex_energy[:125], cld_global[:125], label='RIPL Parameterization')
    plt.plot(ex_energy[:125], cld_tot[:125], label='Adjusted Parameterization')
    plt.step(x[:125], y[:125], label='Empirical Data')

    plt.grid()
    plt.legend()
    #plt.yscale('log')
    plt.xlabel('Excitation Energy (MeV)')
    plt.ylabel('Cumulative Level Distribition')
    plt.title('Comparison of Fitted Parameters to Empirical Levels of $^{52}$Cr')
    plt.savefig('/home/agolas/Desktop/cld_ct.png', dpi=900)


    #     # print(jpi)
    #print(A_1, A)
    #     print(B_calc,B)
    #     print(C_calc,C)
    #
    # Bn = 8.635
    # D_0 = 0.55/1000.
    # D_0_inv = 1./D_0
    # D_0_err = 0.1/1000.
    #
    #ex_energy = np.linspace(1., Bn+6)
    # #ex_energy = x
    #rho_Bn_global = (B)*np.exp(B*(x-C))
    # rho_Bn_adj = B_calc*np.exp(B_calc*(ex_energy -C_calc))
    #
    # cld_global = np.exp(B*(ex_energy-C))
    # cld_adj = np.exp(B_calc*(ex_energy -C_calc))
    #
    # plt.plot(ex_energy, cld_global, label='RIPL Parameterization')
    # plt.plot(ex_energy, cld_adj, label='Adjusted Parameterization')
    # plt.step(x, y, label='Empirical Data')
    # plt.grid()
    # plt.legend()
    # #plt.yscale('log')
    # plt.xlabel('Excitation Energy')
    # plt.ylabel('Cumulative Level Distribition')
    # plt.title('Comparison of RIPL Parameters to Fitted Parameters with $^{52}$Cr')
    # plt.savefig('/home/agolas/Desktop/cld_ct.png', dpi=900)
    # # plt.close()
    #
    #
    #
    # plt.errorbar(Bn, D_0_inv, xerr=D_0_err, fmt='*', color='b', label='$D_0^{-1}$')
    # plt.plot(ex_energy, rho_Bn_global, label='RIPL parameterization')
    # plt.plot(ex_energy, rho_Bn_adj, label='Adjusted Parameterization')
    # plt.grid()
    # plt.legend()
    # plt.yscale('log')
    # plt.title('Level Density at Neutron Separation Energy of $^{91}$Zr')
    # plt.xlabel('Excitation Energy (MeV)')
    # plt.ylabel('Level Density (1/MeV)')
    # plt.savefig('/home/agolas/Desktop/D_0_91.png')
    # plt.close()


    #     cld_array = mc.smoothing_2d(window=5, mode='exponential')
    #     cld_est = cld_array[0]
    #     param_array = cld_array[1]
    #
    #     A_est = param_array[:,0]
    #     B_est = param_array[:,1]
    #     C_est = param_array[:,2]
    #
    #     spin_dep_calc = A_est + globalparams[0]
    #     temp_calc = 1/(B_est + globalparams[1])
    #     cutoff_calc = C_est + globalparams[2]
    #     #print(cutoff_calc)
    #
    #     label = jpi
    #     #cld_data[label] = ripl_lda.cld_hist
    #     #cld_estimate[label] = cld_est
    #     #param_estimate[label] = [spin_dep_calc, temp_calc, cutoff_calc]
    #
    #     with open('/home/agolas/data.dat', 'w') as f:
    #         for line in cld_est:
    #             str_line = str(line)
    #             f.write(str_line)
    #             f.write('\n')
    #
    # fig, ax = plt.subplots(nrows=1, ncols=2,
    #         figsize=(12,5))
    # for n, key in enumerate(cld_estimate.keys()):
    #     #x_data, y_data = cld_data[key][:, :150]
    #     x_est = cld_est[0, :]
    #     y_est = cld_est[1, :]
    #
    #     spin_dep, temp, cutoff = param_estimate[key]
    #
    #     ax[0].step(x_data, y_data, label='actual {}'.format(key))
    #     ax[0].semilogy(x_est, y_est, label='Calculated w/ Energy-Dependent Parameters')
    #
    #     ax[1].plot(x_data, spin_dep, label='Spin-Dependence')
    #     ax[1].plot(x_data, temp, label='Temperature')
    #     ax[1].plot(x_data, cutoff, label='Cutoff Energy')
    #
    #     ax[0].set_xlabel('Excitation Energy (MeV)')
    #     ax[1].set_xlabel('Excitation Energy (MeV)')
    #
    #     ax[0].set_ylabel('Cumulative Level Distribution')
    #     ax[1].set_ylabel('Parameter Value (MeV)')
    #
    #     ax[0].legend()
    #     ax[1].legend()
    #
    #     ax[0].grid()
    #     ax[1].grid()
    # #    if key is not 'total':
    # #        jeepee = 'J' + str(key[0]) + 'pi' + str(key[1])
    # #    else:
    # #        jeepee = key
    # fig.savefig('/home/agolas/Pictures/ctfspins.png')
    # plt.close()

def ripl_bsfgm_comparison(target):

    inputdict = {'target' : target}

    try:
        ripl_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='RIPL')
    except:
        return("Target {} does not exist".format(target))

    bsfgm_lda = level_parameters.LevelDensityAnalyzer(inputdict, source='HFB')

    A = ripl_lda.mass_number
    Z = ripl_lda.num_protons
    compound_label = ripl_lda.compound_label
    N = A - Z

    print("Attempting to determine RIPL/BSFGM deviation for nucleus {}".format(compound_label))

    with open('level_param.json', 'r') as f:
        data = json.load(f)

    try:
        umax = data[compound_label]["Umax"]
    except:
        umax = 14.0

    if umax == 0.0:
        umax = 14.0

    print("Using UMax = {} for Nucleus {}".format(umax, compound_label))



    jpis = []
    for Js in ripl_lda.Jpi:
        if Js[0] == -1: continue
        if Js[1] == 0: continue
        if Js == 'total': continue
        jpis.append(Js)

    hfb_vals = []
    ripl_vals = []

    try:
        for jpi in jpis:
            rval = ripl_lda.cld_extract(ex_energy=umax, j_pi=jpi)
            hval = hfb_lda.cld_extract(ex_energy=umax, j_pi=jpi)

            ripl_vals.append(rval)
            hfb_vals.append(hval)

        hfb_norm = Math(array=hfb_vals).normalize()
        ripl_norm = Math(array=ripl_vals).normalize()

        jxpis = []

        error_per_jpi = np.abs(hfb_norm - ripl_norm)
        maximum_error = np.max(error_per_jpi)
        cumulative_error = np.sum(error_per_jpi)
    except:
        print("Error was encountered with {}".format(compound_label))
        return(Z, N, A, np.nan, np.nan)
    for js in jpis:
        jxpi = js[0]*js[1]
        jxpis.append(jxpi)


    print('RIPL/HFB Deviation for Compound {}: Max={}, Cumulative={}'.format(compound_label, maximum_error, cumulative_error))

    return (Z, N, A, cumulative_error, maximum_error)




def histogram(target):

    inputdict = {'target' : target}

    try:
        ripl_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='RIPL')
        hfb_lda = experimental_fit.LevelDensityAnalyzer(inputdict, source='HFB')
    except:
        return("Target {} does not exist".format(target))

    A = ripl_lda.mass_number
    Z = ripl_lda.num_protons
    compound_label = ripl_lda.compound_label
    N = A - Z

    print("Attempting to determine RIPL/HFB deviation for nucleus {}".format(compound_label))

    with open('level_param.json', 'r') as f:
        data = json.load(f)

    try:
        umax = data[compound_label]["Umax"]
    except:
        umax = 14.0

    if umax == 0.0:
        umax = 14.0

    print("Using UMax = {} for Nucleus {}".format(umax, compound_label))



    jpis = []
    for Js in ripl_lda.Jpi:
        if Js[0] == -1: continue
        if Js[1] == 0: continue
        if Js == 'total': continue
        jpis.append(Js)

    hfb_vals = []
    ripl_vals = []

    try:
        for jpi in jpis:
            rval = ripl_lda.cld_extract(ex_energy=umax, j_pi=jpi)
            hval = hfb_lda.cld_extract(ex_energy=umax, j_pi=jpi)

            ripl_vals.append(rval)
            hfb_vals.append(hval)

        hfb_norm = Math(array=hfb_vals).normalize()
        ripl_norm = Math(array=ripl_vals).normalize()

        jxpis = []

        error_per_jpi = np.abs(hfb_norm - ripl_norm)
        maximum_error = np.max(error_per_jpi)
        cumulative_error = np.sum(error_per_jpi)
    except:
        print("Error was encountered with {}".format(compound_label))
        return(Z, N, A, np.nan, np.nan)
    for js in jpis:
        jxpi = js[0]*js[1]
        jxpis.append(jxpi)


    title_form = "Normalized Distribution of {}".format(compound_label)


    plt.bar(jxpis, hfb_norm, align='edge', width=-0.4, label='HFB')
    plt.bar(jxpis, ripl_norm, align='edge', width=0.4, label='RIPL')


    plt.xlabel('$J*\pi$')
    plt.ylabel('$CLD/NCUMUL$')
    plt.title(title_form)
    plt.legend()

    plt.grid()
    plt.savefig('/home/agolas/Pictures/{}U14.png'.format(compound_label))
    plt.close()

    print('RIPL/HFB Deviation for Compound {}: Max={}, Cumulative={}'.format(compound_label, maximum_error, cumulative_error))

    return (Z, N, A, cumulative_error, maximum_error)


def knownvsunknown(target):

    inputdict = {'target' : target}

    pldl = experimental_fit.PostLevelDensityLoad(inputdict, source='RIPL')
    Jpi = pldl.Jpi


    jpis = []
    for Js in Jpi:
        if Js[0] == -1: continue
        if Js[1] == 0: continue
        if Js == 'total': continue
        jpis.append(Js)

    freq = len(jpis)





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

def main():

    data = np.loadtxt('/home/agolas/Downloads/ZNupdated.dat')

    ZN_grid = np.zeros((71,55))
    ZN_grid[:] = np.nan
    A_list = []
    cumulative_list = []
    maximum_list = []

    magic_numbers = np.asarray([2, 8, 20, 28, 50, 82, 126])
    Z_magic = np.asarray([2, 6, 14, 28, 50, 82, 114])
    N_magic = np.asarray([2, 8, 14, 28, 50, 82, 126])

    dm = []
    dm_0 = []
    dm_1 = []
    dm_2 = []
    dm_3 = []
    dm_4 = []
    dm_5 = []
    dm_6 = []
    dm_7 = []
    dm_8 = []
    dm_9 = []
    dm_10 = []
    dm_11 = []
    dm_12 = []
    dm_13 = []
    dm_14 = []
    nm = []

    A_dm = []
    A_dm_0 = []
    A_dm_1 = []
    A_dm_2 = []
    A_dm_3 = []
    A_dm_4 = []
    A_dm_5 = []
    A_dm_6 = []
    A_dm_7 = []
    A_dm_8 = []
    A_dm_9 = []
    A_dm_10 = []
    A_dm_11 = []
    A_dm_12 = []
    A_dm_13 = []
    A_dm_14 = []
    A_nm = []


    ZN_max_err = ZN_grid
    ZN_cum_err = ZN_grid
    ZN_dist = ZN_grid

    npf = []
    dev = []

    for line in data:
        x = int(line[0])
        y = int(line[1])
        A = int(line[2])
        ZN_max_err[y,x] = line[3]
        ZN_cum_err[y,x] = line[4]

        if line[3] == np.nan: continue
        if line[3] == 0: continue
        #A_list.append(int(line[2]))
        #cumulative_list.append(line[4])
        #maximum_list.append(line[3])

        Ztest = np.abs(magic_numbers - x)
        Ntest = np.abs(magic_numbers - y)
        #
        #
        Zdiff = float(np.min(Ztest))
        Ndiff = float(np.min(Ntest))
        #
        #
        #
        dm_distance = np.sqrt(Zdiff**2 + Ndiff**2)
        #print(np.min(Ztest), np.min(Ntest), int(line[2]))
        #print(x, y, dm_distance)
        

        Zshell = magic_numbers - x
        Nshell = magic_numbers - y

        Zshell = Zshell[Zshell>=0]
        Nshell = Nshell[Nshell>=0]

        N_p = np.min(Zshell)
        N_n = np.min(Nshell)

        if N_n + N_p == 0:
            nucleonic_promiscuity_factor = -1
        else:
            nucleonic_promiscuity_factor =  #(N_n*N_p)/(N_n+N_p)

        ZN_dist[y,x] = -nucleonic_promiscuity_factor

        if nucleonic_promiscuity_factor == -1:
            A_dm.append(int(line[2]))
            dm.append(line[3])

            #else:
            #    Amagic_Z.append(int(line[2]))
            #    magic_Z.append(line[3])

        elif nucleonic_promiscuity_factor == 0:
            A_dm_0.append(int(line[2]))
            dm_0.append(line[3])

        elif nucleonic_promiscuity_factor <= 1:
            A_dm_1.append(int(line[2]))
            dm_1.append(line[3])


        elif nucleonic_promiscuity_factor <=2:
            A_dm_2.append(int(line[2]))
            dm_2.append(line[3])


        elif nucleonic_promiscuity_factor <=3:
            A_dm_3.append(int(line[2]))
            dm_3.append(line[3])


        elif nucleonic_promiscuity_factor <=4:
            A_dm_4.append(int(line[2]))
            dm_4.append(line[3])


        elif nucleonic_promiscuity_factor <=5:
            A_dm_5.append(int(line[2]))
            dm_5.append(line[3])


        elif nucleonic_promiscuity_factor <=6:
            A_dm_6.append(int(line[2]))
            dm_6.append(line[3])


        elif nucleonic_promiscuity_factor <=7:
            A_dm_7.append(int(line[2]))
            dm_7.append(line[3])


        elif nucleonic_promiscuity_factor <=8:
            A_dm_8.append(int(line[2]))
            dm_8.append(line[3])


        elif nucleonic_promiscuity_factor <=9:
            A_dm_9.append(int(line[2]))
            dm_9.append(line[3])

        elif nucleonic_promiscuity_factor <=10:
            A_dm_10.append(int(line[2]))
            dm_10.append(line[3])

        elif nucleonic_promiscuity_factor <=11:
            A_dm_11.append(int(line[2]))
            dm_11.append(line[3])

        elif nucleonic_promiscuity_factor <=12:
            A_dm_12.append(int(line[2]))
            dm_12.append(line[3])


        elif nucleonic_promiscuity_factor <=13:
            A_dm_13.append(int(line[2]))
            dm_13.append(line[3])

        elif nucleonic_promiscuity_factor <=14:
            A_dm_14.append(int(line[2]))
            dm_14.append(line[3])

        npf.append(nucleonic_promiscuity_factor)
        dev.append(line[3])




    dm = np.asarray(dm)
    dm_dev_ave = np.nanmean(dm)
    dm_dev_var = np.nanvar(dm)

    dm_0 = np.asarray(dm_0)
    dm_0_dev_ave = np.nanmean(dm_0)
    dm_0_dev_var = np.nanvar(dm_0)

    dm_1 = np.asarray(dm_1)
    dm_1_dev_ave = np.nanmean(dm_1)
    dm_1_dev_var = np.nanvar(dm_1)

    dm_2 = np.asarray(dm_2)
    dm_2_dev_ave = np.nanmean(dm_2)
    dm_2_dev_var = np.nanvar(dm_2)

    dm_3 = np.asarray(dm_3)
    dm_3_dev_ave = np.nanmean(dm_3)
    dm_3_dev_var = np.nanvar(dm_3)

    dm_4 = np.asarray(dm_4)
    dm_4_dev_ave = np.nanmean(dm_4)
    dm_4_dev_var = np.nanvar(dm_4)

    dm_5 = np.asarray(dm_5)
    dm_5_dev_ave = np.nanmean(dm_5)
    dm_5_dev_var = np.nanvar(dm_5)

    dm_6 = np.asarray(dm_6)
    dm_6_dev_ave = np.nanmean(dm_6)
    dm_6_dev_var = np.nanvar(dm_6)

    dm_7 = np.asarray(dm_7)
    dm_7_dev_ave = np.nanmean(dm_7)
    dm_7_dev_var = np.nanvar(dm_7)

    dm_8 = np.asarray(dm_8)
    dm_8_dev_ave = np.nanmean(dm_8)
    dm_8_dev_var = np.nanvar(dm_8)

    dm_9 = np.asarray(dm_9)
    dm_9_dev_ave = np.nanmean(dm_9)
    dm_9_dev_var = np.nanvar(dm_9)

    dm_10 = np.asarray(dm_10)
    dm_10_dev_ave = np.nanmean(dm_10)
    dm_10_dev_var = np.nanvar(dm_10)

    dm_11 = np.asarray(dm_11)
    dm_11_dev_ave = np.nanmean(dm_11)
    dm_11_dev_var = np.nanvar(dm_11)

    dm_12 = np.asarray(dm_12)
    dm_12_dev_ave = np.nanmean(dm_12)
    dm_12_dev_var = np.nanvar(dm_12)

    dm_13 = np.asarray(dm_13)
    dm_13_dev_ave = np.nanmean(dm_13)
    dm_13_dev_var = np.nanvar(dm_13)

    dm_14 = np.asarray(dm_14)
    dm_14_dev_ave = np.nanmean(dm_14)
    dm_14_dev_var = np.nanvar(dm_14)

    #dm_13 = sorted(dm_13)
    #dm_13 = np.asarray(dm_13[:-1])

    data = [dm, dm_0, dm_1, dm_2, dm_3, dm_4, dm_5, dm_6, dm_7, dm_8,
                    dm_9, dm_10, dm_11, dm_12, dm_13, dm_14]

    means = [dm_dev_ave, dm_0_dev_ave, dm_1_dev_ave, dm_2_dev_ave, dm_3_dev_ave,
            dm_4_dev_ave, dm_5_dev_ave, dm_6_dev_ave, dm_7_dev_ave, dm_8_dev_ave,
            dm_9_dev_ave, dm_10_dev_ave, dm_11_dev_ave, dm_12_dev_ave,
            dm_13_dev_ave, dm_14_dev_ave]
    data_to_plot = []

    for devs in data:
        devs = devs[~np.isnan(devs)]
        data_to_plot.append(devs)


    plt.boxplot(data_to_plot, positions=range(-1,15))


    # plt.errorbar(14, dm_14_dev_ave, yerr=dm_14_dev_var, fmt='o', color='b')
    # plt.errorbar(13, dm_13_dev_ave, yerr=dm_13_dev_var, fmt='o', color='b')
    # plt.errorbar(12, dm_12_dev_ave, yerr=dm_12_dev_var, fmt='o', color='b')
    # plt.errorbar(11, dm_11_dev_ave, yerr=dm_11_dev_var, fmt='o', color='b')
    # plt.errorbar(10, dm_10_dev_ave, yerr=dm_10_dev_var, fmt='o', color='b')
    # plt.errorbar(9, dm_9_dev_ave, yerr=dm_9_dev_var, fmt='o', color='b')
    # plt.errorbar(8, dm_8_dev_ave, yerr=dm_8_dev_var, fmt='o', color='b')
    # plt.errorbar(7, dm_7_dev_ave, yerr=dm_7_dev_var, fmt='o', color='b')
    #
    # plt.errorbar(6, dm_6_dev_ave, yerr=dm_6_dev_var, fmt='o', color='b')
    # plt.errorbar(5, dm_5_dev_ave, yerr=dm_5_dev_var, fmt='o', color='b')
    # plt.errorbar(4, dm_4_dev_ave, yerr=dm_4_dev_var, fmt='o', color='b')
    # plt.errorbar(3, dm_3_dev_ave, yerr=dm_3_dev_var, fmt='o', color='b')
    # plt.errorbar(2, dm_2_dev_ave, yerr=dm_2_dev_var, fmt='o', color='b')
    # plt.errorbar(1, dm_1_dev_ave, yerr=dm_1_dev_var, fmt='o', color='b')
    # plt.errorbar(0, dm_0_dev_ave, yerr=dm_0_dev_var, fmt='o', color='b')
    # plt.errorbar(-1, dm_dev_ave, yerr=dm_dev_var, fmt='o', color='b')

    plt.ylabel('Deviation')
    plt.xlabel('Nucleonic Promiscuity Factor')
    plt.grid(axis='y')
    #plt.xticks(range(1,16), range(-1,15))
    plt.title('Cumulative Deviation between RIPL and HFB')
    plt.savefig('/home/agolas/Pictures/A_dist_error_max_dev_A.png', dpi=900)
    plt.close()


    # #plt.scatter(A_nm, nm, label="DM-(9+)")
    # plt.scatter(A_dm_9, dm_9, label="DM-9")
    # plt.scatter(A_dm_8, dm_8, label="DM-8")
    # plt.scatter(A_dm_7, dm_7, label="DM-7")
    # plt.scatter(A_dm_6, dm_6, label="DM-6")
    # plt.scatter(A_dm_5, dm_5, label="DM-5")
    # plt.scatter(A_dm_4, dm_4, label="DM-4")
    # plt.scatter(A_dm_3, dm_3, label="DM-3")
    # plt.scatter(A_dm_2, dm_2, label="DM-2")
    # plt.scatter(A_dm_1, dm_1, label="DM-1")
    # plt.scatter(A_dm, dm, label="DM")
    #

    plt.scatter(npf, dev)
    plt.title('Maximum Deviation between RIPL and HFB as a Function of P')
    plt.xlabel('Nucleonic Promiscuity Factor')
    plt.ylabel('Maximum Deviation')
    plt.grid()
    plt.legend()
    plt.savefig('/home/agolas/Pictures/sorted_max_dev_A.png', dpi=900)
    plt.close()

    plt.pcolormesh(ZN_dist, vmin=None, vmax=None, cmap='jet')#norm=matplotlib.colors.LogNorm())
    plt.title('Nucleonic Promiscuity Factor')
    plt.ylabel('N')
    plt.xlabel('Z')
    plt.colorbar()
    plt.grid()
    plt.savefig('/home/agolas/Pictures/nucleus_with_A_dist.png', dpi=900)
    plt.close()




    # plt.pcolormesh(ZN_max_err, vmin=None, vmax=None, cmap='jet')#norm=matplotlib.colors.LogNorm())
    # plt.title('Maximum Deviation between RIPL and HFB')
    # plt.ylabel('N')
    # plt.xlabel('Z')
    # plt.colorbar()
    # plt.grid()
    # plt.savefig('/home/agolas/Pictures/u14_max_dev.png', dpi=900)
    # plt.close()
    #
    # plt.pcolormesh(ZN_cum_err, vmin=None, vmax=None, cmap='jet')
    # plt.title('Cumulative Deviation between RIPL and HFB')
    # plt.ylabel('N')
    # plt.xlabel('Z')
    # plt.colorbar()
    # plt.grid()
    # plt.savefig('/home/agolas/Pictures/u14_cum_dev.png', dpi=900)
    # plt.close()
    #
    # plt.scatter(A_list, cumulative_list)
    # plt.title('Cumulative Deviation between RIPL and HFB as a Function of A')
    # plt.xlabel('A')
    # plt.ylabel('Deviation')
    # plt.grid()
    # plt.savefig('/home/agolas/Pictures/u14_cum_dev_A.png', dpi=900)
    # plt.close()
    #
    # plt.scatter(A_list, maximum_list)
    # plt.title('Maximum Deviation between RIPL and HFB as a Function of A')
    # plt.xlabel('A')
    # plt.ylabel('Deviation')
    # plt.grid()
    # plt.savefig('/home/agolas/Pictures/u14_max_dev_A.png', dpi=900)
    # plt.close()




if __name__ == '__main__':
    run()
#
#     Z_label = {
#                1: 'H',
#                2: 'He',
#                3: 'Li',
#                4: 'Be',
#                5: 'B',
#                6: 'C',
#                7: 'N',
#                8: 'O',
#                9: 'F',
#                10: 'Ne',
#                11: 'Na',
#                12: 'Mg',
#                13: 'Al',
#                14: 'Si',
#                15: 'P',
#                16: 'S',
#                17: 'Cl',
#                18: 'Ar',
#                19: 'K',
#                20: 'Ca',
#                21: 'Sc',
#                22: 'Ti',
#                23: 'V',
#                24: 'Cr',
#                25: 'Mn',
#                26: 'Fe',
#                27: 'Co',
#                28: 'Ni',
#                29: 'Cu',
#                30: 'Zn',
#                31: 'Ga',
#                32: 'Ge',
#                33: 'As',
#                34: 'Se',
#                35: 'Br',
#                36: 'Kr',
#                37: 'Rb',
#                38: 'Sr',
#                39: 'Y',
#                40: 'Zr',
#                41: 'Nb',
#                42: 'Mo',
#                43: 'Tc',
#                44: 'Ru',
#                45: 'Rh',
#                46: 'Pd',
#                47: 'Ag',
#                48: 'Cd',
#                49: 'In',
#                50: 'Sn',
#                51: 'Sb',
#                52: 'Te',
#                53: 'I',
#                54: 'Xe'
#                 }

#     #
#     # A_x = []
#     target_array = []
#     # umax_dict = {}
#     #
#     #with open('level_param.json', 'r') as f:
#     #     data = json.load(f)
#     #
#     data = label.split('\n')
#     for line in data:
#         sline = line.split()
#         a_label = str(int(sline[2])-1)
#         Z = int(sline[0])
#         name = Z_label[Z]
#         target = a_label + name
#         target_array.append(target)
#
#     pool = mp.Pool(3)
#     out = pool.map(histogram, target_array)
#     for line in out:
#         val = '{}  {}  {}  {}  {}'.format(line[0], line[1], line[2], line[3], line[4])
#         print(val)
    #np.savetxt('/home/agolas/Documents/ZNOxygen.dat', out)

    #results_list = []

    # for nucleus in target_array:
    #     try:
    #         result = pool.apply_async(histogram, args=(nucleus,))
    #         print(result.get())
    #         results_list.append(result.get())
    #     except:
    #         print('Compound {} Failed'.format(nucleus))


    # pool = mp.Pool(3)
    # out = pool.map(histogram, target_array)
    # print(out)
    # with open('ZNarray.dat', 'a') as f:
    #     f.write(out)
