import json

with open('levels-param.dat', 'r') as f:
    data = list(f)

datafile = {}
for line in data:
    if line.startswith('#'): continue
    sline = line.split()
    Z = None
    A = None
    T = None
    dT = None
    U0 = None
    dU0 = None
    Nlev = None
    Nmax = None
    N0 = None
    Nc = None
    Umax = None
    Uc = None
    Chi = None
    sigma = None

    key = sline[1] + sline[2]
    Z = int(sline[0])
    A = int(sline[1])
    T = float(sline[3])
    dT = float(sline[4])
    U0 = float(sline[5])
    dU0 = float(sline[6])
    Nlev = float(sline[7])
    Nmax = float(sline[8])
    N0 = float(sline[9])
    Nc = float(sline[10])
    Umax = float(sline[11])
    Uc = float(sline[12])

    datafile[key] = { 'T' : T,
                     'dT' : T,
                     'U0' : U0,
                     'dU0' : dU0,
                     'Nlev' : Nlev,
                     'Nmax' : Nmax,
                     'N0' : N0,
                     'Nc' : Nc,
                     'Umax' : Umax,
                     'Uc' : Uc}

#print(datafile)
with open('level_param.json', 'w') as f:
        json.dump(datafile, f, indent=4)
