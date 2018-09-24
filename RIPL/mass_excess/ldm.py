import json

with open('mass_mev.json', 'r') as f:
    data = json.load(f)

#shell_correction = {}

for isotope in data:
        mass = data[isotope]['mass']
        A = data[isotope]['A']
        Z = data[isotope]['Z']
        N = A - Z

        k = 1.79
        a = [15.677, 18.56]
        ci = (1-k*((N-Z)/A)**2)
        c = [a[0]*ci, a[1]*ci, 0.717, 1.21129]

        evol = -c[0]*A
        esur = c[1]*A**(0.66667)
        ecoul = c[2]*(Z**2.0/(A**(0.33333))) - c[3]*(Z**2.0/A)

        if A%2 == 0:
            if Z%2 == 0:
                delta_m = -11.0/(A**(0.5))
            else:
                delta_m = 11.0/(A**(0.5))
        else:
            delta_m = 0.0

        m_n = 8.07144
        m_h = 7.28899

        m_ldm = m_n*N + m_h*Z + evol + esur + ecoul + delta_m
        print(m_ldm)
#        delta_w = mass - m_ldm

#        shell_correction[isotope] = delta_w

#with open('delta_w.json', 'w') as f:
#        json.dump(shell_correction, f, indent=4)
