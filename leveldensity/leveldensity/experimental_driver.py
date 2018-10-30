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

def main():

    data = np.loadtxt('ZNarray3.dat')

    ZN_grid = np.zeros((71,55))
    ZN_grid[:] = np.nan

    ZN_max_err = ZN_cum_err = ZN_grid

    for line in data:
        x = int(line[0])
        y = int(line[1])
        ZN_max_err[y,x] = line[3]
        ZN_cum_err[y,x] = line[4]

    plt.pcolormesh(ZN_max_err)
    plt.title('Maximum Deviation between RIPL and HFB')
    plt.ylabel('N')
    plt.xlabel('Z')
    plt.colorbar()
    plt.grid()
    plt.savefig('/home/agolas/Pictures/max_dev.png', dpi=900)
    plt.close()

    plt.pcolormesh(ZN_cum_err)
    plt.title('Cumulative Deviation between RIPL and HFB')
    plt.ylabel('N')
    plt.xlabel('Z')
    plt.colorbar()
    plt.grid()
    plt.savefig('/home/agolas/Pictures/cum_dev.png', dpi=900)
    plt.close()




if __name__ == '__main__':

    Z_label = {
               1: 'H',
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
               54: 'Xe'
                }
    #
    label = """27  20  47   nan   nan
27  21  48   nan   nan
27  22  49   nan   nan
27  23  50   nan   nan
27  24  51   nan   nan
27  25  52   nan   nan
27  26  53   nan   nan
27  27  54   nan   nan
27  28  55  0.70877937093284  0.1606608080800147
27  29  56  0.8016534401433815  0.1195315273361118
27  30  57  0.873081053451585  0.15346166038793008
27  31  58  0.6371540317156764  0.12485773288976544
27  32  59  1.0417244405843409  0.1481040371525419
27  33  60  0.6377286831545174  0.1039053810684218
27  34  61  0.6907702956327285  0.16102312365096857
27  35  62  0.8720298597351159  0.2884421049698712
27  36  63  0.4802138802775077  0.18037668275140226
27  37  64  0.6533629938025947  0.2672285828962496
27  38  65   nan   nan
27  39  66   nan   nan
27  40  67   nan   nan
27  41  68   nan   nan
27  42  69   nan   nan
27  43  70   nan   nan
27  44  71   nan   nan
27  45  72   nan   nan
27  46  73   nan   nan
27  47  74   nan   nan
27  48  75   nan   nan
28  20  48   nan   nan
28  21  49   nan   nan
28  22  50   nan   nan
28  23  51   nan   nan
28  24  52   nan   nan
28  25  53   nan   nan
28  26  54   nan   nan
28  27  55   nan   nan
28  28  56  1.0082645926817477  0.17911233893671108
28  29  57  0.7052990151377745  0.10290872033105683
28  30  58  1.2695455632642998  0.177633609293392
28  31  59  1.1141836459501577  0.14458001657711517
28  32  60  0.799987745858557  0.10544992702643022
28  33  61  0.49514795728601524  0.12372081903868459
28  34  62  0.39489273204450226  0.09617866323641208
28  35  63  0.3806363055109556  0.07398579343332887
28  36  64   nan   nan
28  37  65   nan   nan
28  38  66   nan   nan
28  39  67   nan   nan
28  40  68  0.7999999999999999  0.39999999999999997
28  41  69   nan   nan
28  42  70   nan   nan
28  43  71   nan   nan
28  44  72   nan   nan
28  45  73   nan   nan
28  46  74   nan   nan
28  47  75   nan   nan
28  48  76   nan   nan
28  49  77   nan   nan
28  50  78   nan   nan
29  23  52   nan   nan
29  24  53   nan   nan
29  25  54   nan   nan
29  26  55   nan   nan
29  27  56   nan   nan
29  28  57   nan   nan
29  29  58  1.1725804390197143  0.1820372160270953
29  30  59  0.8826037722113449  0.11381969314568452
29  31  60  0.90835792495122  0.15501244121158358
29  32  61  0.8012145633636568  0.1300387352368548
29  33  62  0.923076923076923  0.2108318384215302
29  34  63  0.47175772293495677  0.09417625436804178
29  35  64  0.7046448607485862  0.12925905661215437
29  36  65  0.5152223985957234  0.1225178716962891
29  37  66  0.5093580514045031  0.11600693775084293
29  38  67  0.3264114361207066  0.09099033207645388
29  39  68   nan   nan
29  40  69  0.38075636971124643  0.10076340863936216
29  41  70   nan   nan
29  42  71   nan   nan
29  43  72   nan   nan
29  44  73   nan   nan
29  45  74   nan   nan
29  46  75   nan   nan
29  47  76   nan   nan
29  48  77   nan   nan
29  49  78   nan   nan
29  50  79   nan   nan
29  51  80   nan   nan
29  52  81   nan   nan
30  24  54   nan   nan
30  25  55   nan   nan
30  26  56   nan   nan
30  27  57   nan   nan
30  28  58   nan   nan
30  29  59   nan   nan
30  30  60  0.8522610659666142  0.1816658846224482
30  31  61  0.8402299843541529  0.22809341434618596
30  32  62  1.3430801871786398  0.18398346033418522
30  33  63  0.6028335700159433  0.10661625181868077
30  34  64  0.7452730234435913  0.1513547091987984
30  35  65  0.8154790605881838  0.12738398049731894
30  36  66  0.6015546817585886  0.13768290535373676
30  37  67  0.5857816514943057  0.08385837101740415
30  38  68  0.526770185689929  0.11842381242666075
30  39  69   nan   nan
30  40  70   nan   nan
30  41  71  0.7449676369780831  0.20224406236605058
30  42  72  0.24719443597772178  0.12359721798886092
30  43  73   nan   nan
30  44  74   nan   nan
30  45  75   nan   nan
30  46  76   nan   nan
30  47  77   nan   nan
30  48  78   nan   nan
30  49  79   nan   nan
30  50  80   nan   nan
30  51  81   nan   nan
30  52  82   nan   nan
30  53  83   nan   nan
30  54  84   nan   nan
31  25  56   nan   nan
31  26  57   nan   nan
31  27  58   nan   nan
31  28  59   nan   nan
31  29  60   nan   nan
31  30  61   nan   nan
31  31  62   nan   nan
31  32  63   nan   nan
31  33  64  0.8267021684991778  0.17077335309511388
31  34  65  1.5384615384615383  0.279025102408187
31  35  66  0.5509392093247389  0.09239446101088919
31  36  67  0.6635333035423778  0.08508916542000179
31  37  68  0.6430219950756001  0.10715886992985421
31  38  69  0.5306370055495865  0.1231457198396411
31  39  70  0.718299928819144  0.21078888100843496
31  40  71  0.8061643420759148  0.19750970709083426
31  41  72  0.704869707893552  0.13017808787543605
31  42  73   nan   nan
31  43  74   nan   nan
31  44  75  0.4815944609287318  0.10088785006169226
31  45  76   nan   nan
31  46  77   nan   nan
31  47  78  1.4285714285714284  0.5925906363316944
31  48  79   nan   nan
31  49  80   nan   nan
31  50  81   nan   nan
31  51  82   nan   nan
31  52  83   nan   nan
31  53  84   nan   nan
31  54  85   nan   nan
31  55  86   nan   nan
32  26  58   nan   nan
32  27  59   nan   nan
32  28  60   nan   nan
32  29  61   nan   nan
32  30  62   nan   nan
32  31  63   nan   nan
32  32  64   nan   nan
32  33  65   nan   nan
32  34  66   nan   nan
32  35  67  0.7622270087610473  0.1347768238202245
32  36  68  0.9670014452281342  0.13454886951830375
32  37  69  0.4087140412594261  0.10307201086964311
32  38  70  0.4459803059008075  0.086864882999669
32  39  71  0.8715487291548398  0.17210458854067984
32  40  72  0.4953897735268108  0.05816118408144529
32  41  73   nan   nan
32  42  74  0.4506065069186136  0.09343421205531777
32  43  75   nan   nan
32  44  76  0.35705572750166154  0.10652189781209505
32  45  77  0.33385127221654515  0.14609920473828536
32  46  78  0.6745256231195025  0.15568611482913786
32  47  79   nan   nan
32  48  80   nan   nan
32  49  81   nan   nan
32  50  82   nan   nan
32  51  83   nan   nan
32  52  84   nan   nan
32  53  85   nan   nan
32  54  86   nan   nan
32  55  87   nan   nan
32  56  88   nan   nan
32  57  89   nan   nan
33  27  60   nan   nan
33  28  61   nan   nan
33  29  62   nan   nan
33  30  63   nan   nan
33  31  64   nan   nan
33  32  65   nan   nan
33  33  66   nan   nan
33  34  67   nan   nan
33  35  68  1.538461538461538  0.536666957052008
33  36  69  0.7577236043121773  0.19666062458896055
33  37  70  0.7418265529358854  0.1653355515409301
33  38  71  0.9180141220657918  0.20246287333264737
33  39  72   nan   nan
33  40  73  0.885777223893464  0.28295470614314433
33  41  74  0.874410869023876  0.20071857035265234
33  42  75   nan   nan
33  43  76  0.4817580949199456  0.14922331718723902
33  44  77  0.6327749978022285  0.08703131976836792
33  45  78  0.6329379156853225  0.1611582079545475
33  46  79  0.41217707870546105  0.14754858837017487
33  47  80  nan  nan
33  48  81  0.8333333333333334  0.25
33  49  82   nan   nan
33  50  83  0.5049719382947367  0.18719854558290133
33  51  84   nan   nan
33  52  85   nan   nan
33  53  86   nan   nan
33  54  87   nan   nan
33  55  88   nan   nan
33  56  89   nan   nan
33  57  90   nan   nan
33  58  91   nan   nan
34  30  64   nan   nan
34  31  65   nan   nan
34  32  66   nan   nan
34  33  67   nan   nan
34  34  68   nan   nan
34  35  69   nan   nan
34  36  70  0.9227743387895784  0.25717127874636125
34  37  71   nan   nan
34  38  72   nan   nan
34  39  73  1.255474631550553  0.2850609382280096
34  40  74  0.6032302730582539  0.07398106645931879
34  41  75  1.008153441450909  0.12512356924225865
34  42  76  0.5488955653328417  0.07036626670298965
34  43  77  0.6647614522232337  0.07409037788005393
34  44  78  0.5529984865144625  0.0913979570790756
34  45  79   nan   nan
34  46  80  0.5757551082189317  0.21598126755790495
34  47  81  0.3101159456338083  0.13671728616103418
34  48  82  0.37911871292821475  0.09645074114453411
34  49  83   nan   nan
34  50  84  0.6674568869347406  0.13767161186490595
34  51  85   nan   nan
34  52  86   nan   nan
34  53  87   nan   nan
34  54  88   nan   nan
34  55  89   nan   nan
34  56  90   nan   nan
34  57  91   nan   nan
34  58  92   nan   nan
34  59  93   nan   nan
34  60  94   nan   nan
35  32  67   nan   nan
35  33  68   nan   nan
35  34  69   nan   nan
35  35  70   nan   nan
35  36  71  1.3877977814148028  0.23730450991089283
35  37  72  1.2538230632095775  0.22437558980255887
35  38  73  1.376834986192713  0.22454672394794714
35  39  74   nan   nan
35  40  75  1.495277500919755  0.39066493663525925
35  41  76   nan   nan
35  42  77  1.052861869817179  0.12458441623836833
35  43  78   nan   nan
35  44  79  0.8408347419117376  0.10924185187417729
35  45  80  0.7128411203907431  0.11190749197824132
35  46  81  0.7127945301977482  0.13107058969146088
35  47  82  0.7534341051352343  0.11068803880686195
35  48  83  0.8280856132774865  0.1309569273900717
35  49  84   nan   nan
35  50  85  0.5097606732818404  0.14761670607700905
35  51  86   nan   nan
35  52  87   nan   nan
35  53  88   nan   nan
35  54  89   nan   nan
35  55  90   nan   nan
35  56  91   nan   nan
35  57  92   nan   nan
35  58  93   nan   nan
35  59  94   nan   nan
35  60  95   nan   nan
35  61  96   nan   nan
35  62  97   nan   nan
36  33  69   nan   nan
36  34  70   nan   nan
36  35  71   nan   nan
36  36  72  1.1999999999999997  0.5215245713613503
36  37  73   nan   nan
36  38  74   nan   nan
36  39  75  1.1639357069585121  0.149606094161498
36  40  76   nan   nan
36  41  77  1.2547783377740678  0.2326079878196985
36  42  78   nan   nan
36  43  79   nan   nan
36  44  80  0.8054804496696683  0.12837624427652916
36  45  81  0.9326575855968086  0.1022877856487953
36  46  82  0.47304858299985325  0.06371577565334213
36  47  83  0.9905207899843065  0.25101438811647536
36  48  84   nan   nan
36  49  85  0.6725093035891743  0.11742569598658967
36  50  86  0.6193600867029931  0.12474514127259206
36  51  87  0.8027253995625097  0.3303004020443609
36  52  88  0.43667967216127934  0.10685740463797655
36  53  89   nan   nan
36  54  90  1.0667095759289782  0.4213737017082757
36  55  91  1.1051927919356292  0.2857142857142857
36  56  92   nan   nan
36  57  93   nan   nan
36  58  94   nan   nan
36  59  95   nan   nan
36  60  96   nan   nan
36  61  97   nan   nan
36  62  98   nan   nan
36  63  99   nan   nan
36  64  100   nan   nan
37  34  71   nan   nan
37  35  72   nan   nan
37  36  73   nan   nan
37  37  74   nan   nan
37  38  75   nan   nan
37  39  76   nan   nan
37  40  77   nan   nan
37  41  78   nan   nan
37  42  79   nan   nan
37  43  80  1.1111111111111112  0.18564115118171784
37  44  81  0.7118451522698476  0.09256346397518245
37  45  82   nan   nan
37  46  83  1.1325514546534157  0.14932616881192284
37  47  84  0.7298258828499679  0.12603896106904142
37  48  85  0.8756410721224659  0.12816605872221393
37  49  86  0.9676634694579148  0.1980126507453716
37  50  87  0.38770596967803905  0.13141755503160496
37  51  88  0.6865006597333566  0.13997002292393307
37  52  89   nan   nan
37  53  90   nan   nan
37  54  91   nan   nan
37  55  92   nan   nan
37  56  93   nan   nan
37  57  94   nan   nan
37  58  95   nan   nan
37  59  96   nan   nan
37  60  97   nan   nan
37  61  98   nan   nan
37  62  99   nan   nan
37  63  100   nan   nan
37  64  101   nan   nan
37  65  102   nan   nan
38  35  73   nan   nan
38  36  74   nan   nan
38  37  75   nan   nan
38  38  76   nan   nan
38  39  77   nan   nan
38  40  78   nan   nan
38  41  79   nan   nan
38  42  80   nan   nan
38  43  81   nan   nan
38  44  82  1.0114247598317883  0.19460709207764626
38  45  83  0.8157872881144782  0.1308374928367952
38  46  84  0.8374520176123165  0.12294817055037588
38  47  85  0.5269737002785249  0.07806614173849256
38  48  86  0.6642722219180318  0.15998997711305768
38  49  87   nan   nan
38  50  88  0.6844850391353757  0.08346011411934093
38  51  89  0.5588476696761824  0.1377749074624497
38  52  90  0.42742903138417876  0.15159839710779366
38  53  91  0.6705445702869722  0.18813682588149605
38  54  92  0.6974357025654463  0.14450593381373544
38  55  93  0.6720373804094376  0.1540403513782666
38  56  94   nan   nan
38  57  95  0.9811921800821859  0.264438929099791
38  58  96  0.5497676969055063  0.17329770856343119
38  59  97  0.7534701085169604  0.1419021059266127
38  60  98   nan   nan
38  61  99  0.3144765084664918  0.08911280705938468
38  62  100   nan   nan
38  63  101   nan   nan
38  64  102   nan   nan
38  65  103   nan   nan
38  66  104   nan   nan
38  67  105   nan   nan
38  68  106   nan   nan
39  37  76   nan   nan
39  38  77   nan   nan
39  39  78   nan   nan
39  40  79   nan   nan
39  41  80  1.1599968428253205  0.154865556978233
39  42  81   nan   nan
39  43  82   nan   nan
39  44  83   nan   nan
39  45  84   nan   nan
39  46  85  1.0386498485430005  0.1707612803502052
39  47  86  0.4739292914010718  0.11854968116017454
39  48  87  0.36489846663903847  0.04755491089929548
39  49  88  0.6824036769258839  0.08116078501447554
39  50  89  0.6448850833708878  0.07683641279855791
39  51  90  0.885548744723288  0.14640580661307445
39  52  91  0.5450588819268095  0.08829771162603323
39  53  92  0.46795964762716685  0.23397982381358345
39  54  93   nan   nan
39  55  94  0.35439089065668383  0.15919067673864115
39  56  95   nan   nan
39  57  96   nan   nan
39  58  97  0.7293969414462157  0.22015801085851483
39  59  98   nan   nan
39  60  99  0.48253584524059934  0.05598711062524134
39  61  100   nan   nan
39  62  101  0.396102671642626  0.06484737467026669
39  63  102   nan   nan
39  64  103   nan   nan
39  65  104   nan   nan
39  66  105   nan   nan
39  67  106   nan   nan
39  68  107   nan   nan
39  69  108   nan   nan
40  38  78   nan   nan
40  39  79   nan   nan
40  40  80   nan   nan
40  41  81   nan   nan
40  42  82   nan   nan
40  43  83   nan   nan
40  44  84  1.3906297610591882  0.2099473561942665
40  45  85   nan   nan
40  46  86   nan   nan
40  47  87   nan   nan
40  48  88  0.5542077371713268  0.10105661829701026
40  49  89   nan   nan
40  50  90  0.688738972496108  0.09403217436429684
40  51  91  0.5077628731222016  0.1141304640800116
40  52  92  0.5130979339216976  0.07125009242933525
40  53  93  0.47280145219835257  0.11292931702331516
40  54  94  0.506081093905999  0.11482653374990034
40  55  95  0.5813378647956722  0.08933122069210786
40  56  96  0.5564009108265708  0.19201653270591854
40  57  97   nan   nan
40  58  98  0.730463679131004  0.1375054691935245
40  59  99  0.6343966077045554  0.13752443987172824
40  60  100   nan   nan
40  61  101  0.505142257944543  0.11155325727208862
40  62  102   nan   nan
40  63  103   nan   nan
40  64  104   nan   nan
40  65  105   nan   nan
40  66  106   nan   nan
40  67  107   nan   nan
40  68  108   nan   nan
40  69  109   nan   nan
40  70  110   nan   nan
41  40  81   nan   nan
41  41  82   nan   nan
41  42  83   nan   nan
41  43  84   nan   nan
41  44  85   nan   nan
41  45  86   nan   nan
41  46  87   nan   nan
41  47  88   nan   nan
41  48  89   nan   nan
41  49  90   nan   nan
41  50  91  0.7056305829356455  0.0863237138351288
41  51  92   nan   nan
41  52  93  0.7009945732503992  0.09072197758485992
41  53  94  0.4200867263301697  0.06548761994390337
41  54  95   nan   nan
41  55  96   nan   nan
41  56  97   nan   nan
41  57  98   nan   nan
41  58  99   nan   nan
41  59  100   nan   nan
41  60  101  0.6174713462678535  0.08222132918204136
41  61  102   nan   nan
41  62  103   nan   nan
41  63  104  0.0  0.0
41  64  105   nan   nan
41  65  106   nan   nan
41  66  107   nan   nan
41  67  108   nan   nan
41  68  109   nan   nan
41  69  110   nan   nan
41  70  111   nan   nan
42  41  83   nan   nan
42  42  84   nan   nan
42  43  85   nan   nan
42  44  86   nan   nan
42  45  87   nan   nan
42  46  88   nan   nan
42  47  89   nan   nan
42  48  90   nan   nan
42  49  91  0.6862260473987134  0.10328378086292044
42  50  92  0.6555881431703637  0.08881374754406447
42  51  93  0.3583572199885221  0.06472474404795885
42  52  94  0.5298959945841276  0.07540338587330968
42  53  95  0.6574850477220989  0.11054597263786389
42  54  96  0.4659375802617263  0.06294312389621176
42  55  97  0.539614114219349  0.11907217754069957
42  56  98  0.3746939744418216  0.10951455286375564
42  57  99  0.5729064494215663  0.12610240368896547
42  58  100  0.5833564532515253  0.07662556311138875
42  59  101  0.5319614181036149  0.14364931168145467
42  60  102  0.5485480260678955  0.10947235762042098
42  61  103   nan   nan
42  62  104  0.8242431061523894  0.1382545134932257
42  63  105   nan   nan
42  64  106   nan   nan
42  65  107   nan   nan
42  66  108   nan   nan
42  67  109   nan   nan
42  68  110   nan   nan
42  69  111   nan   nan
42  70  112   nan   nan
43  42  85   nan   nan
43  43  86   nan   nan
43  44  87   nan   nan
43  45  88   nan   nan
43  46  89   nan   nan
43  47  90   nan   nan
43  48  91   nan   nan
43  49  92   nan   nan
43  50  93   nan   nan
43  51  94   nan   nan
43  52  95  1.0216388953913162  0.1358504820071929
43  53  96  0.7151674370381388  0.14800632819890047
43  54  97  0.4647938815623044  0.05953888441640931
43  55  98   nan   nan
43  56  99  0.44102861569815865  0.07244809905547232
43  57  100  0.8492572204897252  0.15942583495840637
43  58  101   nan   nan
43  59  102  0.1908037192714303  0.09540185963571518
43  60  103   nan   nan
43  61  104   nan   nan
43  62  105  1.1428571428571428  0.2207112837086905
43  63  106  0.0  0.0
43  64  107   nan   nan
43  65  108   nan   nan
43  66  109   nan   nan
43  67  110   nan   nan
43  68  111   nan   nan
43  69  112   nan   nan
43  70  113   nan   nan
44  43  87   nan   nan
44  44  88   nan   nan
44  45  89   nan   nan
44  46  90   nan   nan
44  47  91   nan   nan
44  48  92   nan   nan
44  49  93   nan   nan
44  50  94  1.3011406844106463  0.2237642585551331
44  51  95   nan   nan
44  52  96  0.5715185197607588  0.07257930417156688
44  53  97   nan   nan
44  54  98   nan   nan
44  55  99  0.889747493556989  0.12647429499307283
44  56  100  0.617676506821863  0.1017057973038493
44  57  101  0.6775818666287157  0.10489071158447916
44  58  102   nan   nan
44  59  103  0.6015785133702439  0.22674441434012643
44  60  104  0.4534635489259648  0.08454868179575817
44  61  105  0.44900795992499043  0.09709272591368662
44  62  106  0.5781655686454598  0.09633149054171009
44  63  107  0.9111635887009607  0.1519399529890791
44  64  108   nan   nan
44  65  109  0.6467451807581883  0.1191627116657878
44  66  110  0.8886360188135337  0.2646032468517676
44  67  111   nan   nan
44  68  112   nan   nan
44  69  113   nan   nan
44  70  114   nan   nan
45  44  89   nan   nan
45  45  90   nan   nan
45  46  91   nan   nan
45  47  92   nan   nan
45  48  93  0.5559071729957805  0.19672995780590719
45  49  94   nan   nan
45  50  95   nan   nan
45  51  96   nan   nan
45  52  97   nan   nan
45  53  98   nan   nan
45  54  99  0.6800537793094115  0.17292635012342003
45  55  100  1.098399668897013  0.18303525703332563
45  56  101  1.2699615613842277  0.19431595120537593
45  57  102  1.0922299744467596  0.1577800533715885
45  58  103  0.7610444878395514  0.0745291999413149
45  59  104  0.7914066641236209  0.10856732011141546
45  60  105  0.5672686281530591  0.06228281329962623
45  61  106   nan   nan
45  62  107  0.9096176465266091  0.14716718238242588
45  63  108  0.2912883182464307  0.10852336064454986
45  64  109   nan   nan
45  65  110  0.3269534014262911  0.12066513774734887
45  66  111   nan   nan
45  67  112   nan   nan
45  68  113  0.8803633409644783  0.1791859643200017
45  69  114  0.36526196787392384  0.07797821902883814
45  70  115   nan   nan
46  45  91   nan   nan
46  46  92   nan   nan
46  47  93   nan   nan
46  48  94   nan   nan
46  49  95   nan   nan
46  50  96   nan   nan
46  51  97  0.4867163017869482  0.11958187285253838
46  52  98   nan   nan
46  53  99   nan   nan
46  54  100   nan   nan
46  55  101  0.4909309302095126  0.11940536498794338
46  56  102   nan   nan
46  57  103  0.6708271707369864  0.10140326828036506
46  58  104  0.5714367252747554  0.14655798810339196
46  59  105  0.6425520774898088  0.059590859792812007
46  60  106  0.7124544598092638  0.09768469836077912
46  61  107  0.8103277051052983  0.12433382636907561
46  62  108  0.5826923605327512  0.07164481768990183
46  63  109  0.6430904657715312  0.11958895725914151
46  64  110  0.739193783389995  0.15486019565669873
46  65  111  0.7855125134186217  0.1999583969247137
46  66  112   nan   nan
46  67  113   nan   nan
46  68  114   nan   nan
46  69  115   nan   nan
46  70  116  0.7290475206137026  0.2816746671566488
47  46  93   nan   nan
47  47  94   nan   nan
47  48  95   nan   nan
47  49  96   nan   nan
47  50  97   nan   nan
47  51  98  0.6058863427555443  0.2740327086492652
47  52  99   nan   nan
47  53  100   nan   nan
47  54  101   nan   nan
47  55  102   nan   nan
47  56  103  1.05397349750366  0.1403160779308532
47  57  104   nan   nan
47  58  105  1.0553272073783573  0.11545482301594355
47  59  106   nan   nan
47  60  107  1.1461388731548756  0.10309604630282382
47  61  108  0.6139893005831345  0.09618957873405032
47  62  109  0.6732309685466715  0.11300040826001802
47  63  110  0.831585734344632  0.2909932703130167
47  64  111  0.5509980549840465  0.1003432768328395
47  65  112   nan   nan
47  66  113   nan   nan
47  67  114   nan   nan
47  68  115  0.9351329346232387  0.23747797661798078
47  69  116   nan   nan
47  70  117   nan   nan
48  47  95   nan   nan
48  48  96   nan   nan
48  49  97   nan   nan
48  50  98   nan   nan
48  51  99   nan   nan
48  52  100   nan   nan
48  53  101   nan   nan
48  54  102   nan   nan
48  55  103   nan   nan
48  56  104   nan   nan
48  57  105  0.5060975609756099  0.14024390243902446
48  58  106  0.7720158704111475  0.17001413275632
48  59  107  0.725155239493958  0.15632166138465903
48  60  108  0.7748512439484115  0.07425485434202682
48  61  109  0.8678261453908992  0.20817620302180398
48  62  110  0.9272349157610487  0.17418614233161792
48  63  111  0.548578213204712  0.13957710711837032
48  64  112  0.5090206455441144  0.17223255311507144
48  65  113  0.6543980109017883  0.09777692773580696
48  66  114  0.70901376929107  0.18345271506164562
48  67  115  0.6205348763233234  0.10176351978379913
48  68  116  0.4781313284286039  0.06516178045217888
48  69  117  0.7967851242659207  0.2522821276544431
48  70  118  0.2755467196819086  0.10478793903247183
49  48  97   nan   nan
49  49  98   nan   nan
49  50  99   nan   nan
49  51  100   nan   nan
49  52  101   nan   nan
49  53  102   nan   nan
49  54  103   nan   nan
49  55  104   nan   nan
49  56  105   nan   nan
49  57  106   nan   nan
49  58  107  1.161107344387649  0.19181866981884452
49  59  108   nan   nan
49  60  109  1.1932889314957478  0.18167288099239876
49  61  110   nan   nan
49  62  111   nan   nan
49  63  112  0.7561304078645119  0.1333027391553403
49  64  113  1.0133229636160366  0.13420981469158752
49  65  114   nan   nan
49  66  115  0.7651704899577823  0.14423452818195098
49  67  116   nan   nan
49  68  117  0.781276433653651  0.1638350979672775
49  69  118   nan   nan
49  70  119  0.7109333919811736  0.11601961538197833
50  49  99   nan   nan
50  50  100   nan   nan
50  51  101  0.0  0.0
50  52  102   nan   nan
50  53  103   nan   nan
50  54  104   nan   nan
50  55  105   nan   nan
50  56  106   nan   nan
50  57  107  0.8571428571428572  0.23415560884365372
50  58  108  1.3437180876813806  0.26961035238841335
50  59  109  1.0575241807223945  0.34236593466336424
50  60  110  0.5651281123240861  0.13944728124062555
50  61  111  0.9600268740751732  0.19595611117293796
50  62  112  0.6944672631220763  0.1570553823627471
50  63  113  0.8073013016294622  0.14853904460940154
50  64  114  0.7674358593156001  0.10259135897965393
50  65  115  0.7324503866359059  0.05866792891150449
50  66  116  0.6528080531527024  0.07185006358247831
50  67  117  0.9732469610589235  0.1954506699330217
50  68  118   nan   nan
50  69  119  0.7123267566654381  0.13145769141135613
50  70  120  0.45266498625050755  0.06857340194620173
51  52  103   nan   nan
51  53  104   nan   nan
51  54  105   nan   nan
51  55  106   nan   nan
51  56  107   nan   nan
51  57  108   nan   nan
51  58  109  1.2871161214460034  0.2519384509321799
51  59  110   nan   nan
51  60  111   nan   nan
51  61  112  0.6898027523346495  0.19326779302873207
51  62  113  1.284494586725486  0.15256511618201463
51  63  114  1.087075908253029  0.10396647499020338
51  64  115  0.9345664793024431  0.14693353662900055
51  65  116  0.45824810426365303  0.08450434422120534
51  66  117  0.9736345560233965  0.1308544209154497
51  67  118  0.3942108929515528  0.1226820389359219
51  68  119  0.6888390056200966  0.05995623227868319
51  69  120   nan   nan
51  70  121  0.6739924717702408  0.12279893716116237
52  53  105   nan   nan
52  54  106   nan   nan
52  55  107   nan   nan
52  56  108   nan   nan
52  57  109   nan   nan
52  58  110   nan   nan
52  59  111   nan   nan
52  60  112   nan   nan
52  61  113   nan   nan
52  62  114  1.0122595134735521  0.17393775993216798
52  63  115   nan   nan
52  64  116   nan   nan
52  65  117   nan   nan
52  66  118  0.7262516939870423  0.17323778937011533
52  67  119  0.5930085161073344  0.07903142178354514
52  68  120  0.3468094794046047  0.08307169058761019
52  69  121  0.7480596614188462  0.09724818436770671
52  70  122  0.7035459306726355  0.11144099958051124
53  54  107   nan   nan
53  55  108   nan   nan
53  56  109   nan   nan
53  57  110   nan   nan
53  58  111   nan   nan
53  59  112   nan   nan
53  60  113   nan   nan
53  61  114   nan   nan
53  62  115   nan   nan
53  63  116  1.369791522467399  0.25694232952185825
53  64  117  1.3809523809523805  0.1357704433198468
53  65  118   nan   nan
53  66  119  1.1112844497585301  0.0932958669513019
53  67  120   nan   nan
53  68  121  0.791323699774404  0.10968386074019879
53  69  122  1.1488441423044775  0.13280899220620573
53  70  123  0.6394599259090837  0.0843371730346416
54  55  109   nan   nan
54  56  110   nan   nan
54  57  111   nan   nan
54  58  112   nan   nan
54  59  113   nan   nan
54  60  114   nan   nan
54  61  115   nan   nan
54  62  116   nan   nan
54  63  117   nan   nan
54  64  118  0.7228865121701863  0.19343756133428242
54  65  119  0.7560903471325336  0.10548279125711077
54  66  120  1.2317269901200778  0.19063449826103215
54  67  121  1.176221697700848  0.13973249409913455
54  68  122  1.5428571428571427  0.2285735703978175
54  69  123  0.6972960953114551  0.06792805060613513
54  70  124   nan   nan"""
    #
    # A_x = []
    target_array = []
    # umax_dict = {}
    #
    #with open('level_param.json', 'r') as f:
    #     data = json.load(f)
    #
    data = label.split('\n')
    for line in data:
        sline = line.split()
        a_label = str(int(sline[2])-1)
        Z = int(sline[0])
        name = Z_label[Z]
        target = a_label + name
        target_array.append(target)

    pool = mp.Pool(3)
    out = pool.map(histogram, target_array)
    for line in out:
        val = '{}  {}  {}  {}  {}'.format(line[0], line[1], line[2], line[3], line[4])
        print(val)
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
