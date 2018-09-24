import json

with open('delta_w.json', 'r') as f:
    shell_data = json.load(f)

with open('mass_mev.json', 'r') as f:
    mass_data = json.load(f)

with open('level_param.json', 'r') as f:
    levels_data = json.load(f)

with open('/home/agolas/empire/work/leveldensity/data.json', 'r') as f:
    res_data = json.load(f)

datafile = {}
keys = []

for level_isotope in levels_data:
    if level_isotope in mass_data:
        if level_isotope in res_data:
            if level_isotope in shell_data:
                keys.append(level_isotope)
for isotope in keys:
    print(isotope)

    Z     = None
    A     = None
    mass  = None
    T     = None
    dT    = None
    U0    = None
    dU0   = None
    Nlev  = None
    Nmax  = None
    N0    = None
    Nc    = None
    Umax  = None
    Uc    = None
    Bn = None
    Io = None

    Z = mass_data[isotope]['Z']
    A = mass_data[isotope]['A']
    mass = mass_data[isotope]['mass']

    Bn = res_data[isotope]['Bn']
    Io = res_data[isotope]['Io']

    T     = levels_data[isotope]['T']
    dT    = levels_data[isotope]['dT']
    U0    = levels_data[isotope]['U0']
    dU0   = levels_data[isotope]['dU0']
    Nlev  = levels_data[isotope]['Nlev']
    Nmax  = levels_data[isotope]['Nmax']
    N0    = levels_data[isotope]['N0']
    Nc    = levels_data[isotope]['Nc']
    Umax  = levels_data[isotope]['Umax']
    Uc    = levels_data[isotope]['Uc']

    del_w = shell_data[isotope]


    datafile[isotope] = {'Z' : Z,
                         'A' : A,
                         'mass' : mass,
                         'T'    : T,
                         'dT'   : dT,
                         'U0'   : U0,
                         'dU0'  : dU0,
                         'Nlev' : Nlev,
                         'Nmax' : Nmax,
                         'N0'   : N0,
                         'Nc'   : Nc,
                         'Umax' : Umax,
                         'Uc'   : Uc,
                         'Bn'   : Bn,
                         'Io'   : Io,
                         'shell_correction' : del_w}

#print(datafile)
with open('ripl3.json', 'w') as f:
        json.dump(datafile, f, indent=4)
