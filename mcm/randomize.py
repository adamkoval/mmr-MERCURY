# Masses
# 0-15Mj for 2:1, 
# 0-4.1Mj for 5:3,
# 0-3.8Mj for 3:2,
# 0-3Mj for 7:5,
# 0-1.4Mj for 4:3.

import numpy as np
import sys
import os
import re
import argparse

import func as fn

parser = argparse.ArgumentParser()
parser.add_argument('--resonance', '-res',
                    dest='resonance',
                    action='store')
parser.add_argument('--process_no', '-pno',
                    dest='process_no',
                    action='store')
args = parser.parse_args()

# Some definitions
res_str = args.resonance
pno = args.process_no
Tnorm_o = max(float(res_str[0]), float(res_str[1]))
Tnorm_i = min(float(res_str[0]), float(res_str[1]))
Mjup = 1.898e27 # kg
Msol = 1.989e30 # kg

# Mass ranges
Sum = Tnorm_o + Tnorm_i
if Sum == 3: # 2:1
    range0 = 0.
    rangef = 15.
elif Sum == 8: # 5:3
    range0 = 0.
    rangef = 4.1
elif Sum == 5: # 3:2
    range0 = 0.
    rangef = 3.8
elif Sum == 12: # 7:5
    range0 = 0.
    rangef = 3.
elif Sum == 7: # 4:3
    range0 = 0.
    rangef = 1.4
else:
    print("Please input a valid resonance ('21', '53', '32', '75' or '43')")
    sys.exit()

# Masses of inner (planet1) and outer (planet2) planets
new_mi = np.random.uniform(range0, rangef) * Mjup/Msol
new_mo = np.random.uniform(range0, rangef) * Mjup/Msol

# Initial positions (must be improved)
new_a0i = 5

res_float = Tnorm_o / Tnorm_i
out_disp_frac = .4
Dres_out = fn.disp(res_float, out_disp_frac, 'out')
new_a0o = (res_float + Dres_out)**(2/3) * new_a0i

# Migration timescale
new_p1i = 1.e15
new_p1o = 1.e6

# Migration ditance (needs changing)
new_p2i = 0.0

inn_disp_frac = .1
Dres_in = fn.disp(res_float, inn_disp_frac, 'in')
fin_a0o = (res_float - Dres_in)**(2/3) * new_a0i
new_p2o = fin_a0o - new_a0o

# Inclination
new_incli = 0.0
new_inclo = 0.0

# Argument of periapsis
new_arpi = 0.0
new_arpo = np.random.uniform(360)

# Mean anomaly
new_Mi = 0.0
new_Mo = np.random.uniform(360)

####################################
# Amending file
with open('mercury_{}/big.in'.format(pno), 'r') as f:
    bigin = f.read()

# Mass
bigin = re.sub('(planet1.*)m=\S*', '\g<1>m={:.3g}'.format(new_mi), bigin)
bigin = re.sub('(planet2.*)m=\S*', '\g<1>m={:.3g}'.format(new_mo), bigin)

# Migration timescale
bigin = re.sub('(planet1.*)p1=\S*', '\g<1>p1={:.3g}'.format(new_p1i), bigin)
bigin = re.sub('(planet2.*)p1=\S*', '\g<1>p1={:.3g}'.format(new_p1o), bigin)

# Migration distance
bigin = re.sub('(planet1.*)p2=\S*', '\g<1>p2={:.3g}'.format(new_p2i), bigin)
bigin = re.sub('(planet2.*)p2=\S*', '\g<1>p2={:.3g}'.format(new_p2o), bigin)

# Initial position
bigin = re.sub('(planet1.*\n\s*)\S*', '\g<1>{:.3f}'.format(new_a0i), bigin)
bigin = re.sub('(planet2.*\n\s*)\S*', '\g<1>{:.3f}'.format(new_a0o), bigin)

# Inclination
bigin = re.sub('(planet1.*\n\s*\S*\s*\S*\s*)\S*', '\g<1>{:.3f}'.format(new_incli), bigin)
bigin = re.sub('(planet2.*\n\s*\S*\s*\S*\s*)\S*', '\g<1>{:.3f}'.format(new_inclo), bigin)

# Argument of periapsis
bigin = re.sub('(planet1.*\n.*\n\s*)\S*', '\g<1>{:.3f}'.format(new_arpi), bigin)
bigin = re.sub('(planet2.*\n.*\n\s*)\S*', '\g<1>{:.3f}'.format(new_arpo), bigin)

# Mean anomaly
bigin = re.sub('(planet1.*\n.*\n\s*\S*\s*\S*\s*)\S*', '\g<1>{:.3f}'.format(new_Mi), bigin)
bigin = re.sub('(planet2.*\n.*\n\s*\S*\s*\S*\s*)\S*', '\g<1>{:.3f}'.format(new_Mo), bigin)

# Writing new file
new_f = open('mercury_{}/big.in'.format(pno), 'w')
new_f.write(bigin)
new_f.close()
