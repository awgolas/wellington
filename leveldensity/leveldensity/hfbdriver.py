from __future__ import print_function
import BNL.restools.level_density  as ldmodule

RIPLPATH="/home/agolas/empire/RIPL"
levelFile=RIPLPATH+"/levels/z024.dat"
HFBTabFile=RIPLPATH+"/densities/total/level-densities-hfb/z024.tab"
HBFCorFile=RIPLPATH+"/densities/total/level-densities-hfb/z024.cor"

cr52hfb=ldmodule.readHFBMLevelDensityTable(open(HFBTabFile).read(), 24, 52)
Js=cr52hfb.get_J_range()

Emax=3.0
print('CLD at EMax=%s MeV'%(str(Emax)))

with open('hfbdata.dat', 'w') as f:
    init = '{} {}'.format(24, 52)
    f.write(init)
    f.write('\n')
    for Pi in (-1,1):
        print('For Pi=', Pi, 'CLD(3.0 MeV)=', cr52hfb.get_CLD(Pi=Pi).evaluate(17.0))
        print("J", "Pi", "CLD(3.0 MeV)")
        for J in Js:
            if Pi == -1 and J == 0: continue
            try:
                line = '{}  {}  {}'.format(J, Pi,
                        cr52hfb.get_CLD().evaluate(17.0))
                f.write(line)
                f.write('\n')
# FUDGE's level densities use Ex in MeV and rho in 1/MeV
            except KeyError:
                continue
