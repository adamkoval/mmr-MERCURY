import numpy as np
import matplotlib.pyplot as plt

import func as fn

res_str = '53'
path = '../completed/{}/planets'.format(res_str)
planet_i = fn.read_planet('{}/10-planet1.aei'.format(path))
planet_o = fn.read_planet('{}/10-planet2.aei'.format(path))

phi1, phi2, t_phi = fn.resvar(planet_i, planet_o, res_str)
#Lambda_i = planet_i['node'] + planet_i['peri'] # Longidute of periapsis
#Lambda_o = planet_o['node'] + planet_o['peri']
#omega_i = Lambda_i + planet_i['M'] # Mean longitude
#omega_o = Lambda_o + planet_o['M']
#p = int(res_str[0])
#q = int(res_str[1])
#phi1 = p*Lambda_o - q*Lambda_i - (p - q)*omega_o
#phi2 = p*Lambda_o - q*Lambda_i - (p - q)*omega_i


fig, ax = plt.subplots(2)
ax[0].plot(t_phi, phi1)
ax[1].plot(t_phi, phi2)
plt.show()
