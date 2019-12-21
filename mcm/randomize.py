# Masses
# 0-12Mj for 2:1, 
# 0-6Mj for 5:3,
# 0.2-3.8Mj for 3:2,
# 0-3.5Mj for 7:5,
# 0-3.5Mj for 4:3.

import numpy as np
import sys
import os
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--resonance', '-res',
                    dest='resonance',
                    action='store')
args = parser.parse_args()

# Some definitions
resonance = args.resonance
higher = max(float(resonance[0]), float(resonance[1]))
lower = min(float(resonance[0]), float(resonance[1]))
Mjup = 1.898e27 # kg
Msol = 1.989e30 # kg

# Mass ranges
Sum = higher + lower
if Sum == 3:
    range0 = 0.
    rangef = 15.
elif Sum == 8:
    range0 = .2
    rangef = 4.1
elif Sum == 5:
    range0 = 3.5
    rangef = 4.45
elif Sum == 12:
    range0 = 0.93
    rangef = 1.39
elif Sum == 7:
    range0 = 2.
    rangef = 8.
else:
    print("Please input a valid resonance ('21', '53', '32', '75' or '43')")
    sys.exit()

# Masses of inner (planet1) and outer (planet2) planets
new_mi = np.random.uniform(range0, rangef) * Mjup / Msol
new_mo = np.random.uniform(range0, rangef) * Mjup / Msol

# Initial positions (must be improved)
a0_range = (5, 6)
new_a0i = np.random.uniform(*a0_range)

factor = higher/lower
new_a0o = factor**(2/3) * new_a0i + .2

# Migration timescale
new_p1i = 1.e15
new_p1o = 1.e6

# Migration ditance (needs changing)
new_p2i = 0.0
new_p2o = -.3

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
with open('mercury/big.in', 'r') as f:
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
new_f = open('mercury/big.in', 'w')
new_f.write(bigin)
new_f.close()
