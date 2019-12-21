import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import funcs as fn

planets_csv_old = pd.read_csv('planets.csv', skiprows=23, header=0, delim_whitespace=False)
planets_csv_new = pd.read_csv('planets_new.csv', skiprows=34, header=0, delim_whitespace=False)
planets_csv_control = pd.read_csv('planets_control.csv', skiprows=23, header=0, delim_whitespace=False)


resonance = '75'
tol = 15
##observed_old = fn.find_resonances(planets_csv_old, tol)
##observed_new = fn.find_resonances(planets_csv_new, tol)
observed_control, massj, massk = fn.find_resonances(planets_csv_control, tol)

fig, ax, boundary = fn.stability_fig_setup(resonance)
##fn.plot_observed(observed_old, resonance, fig, ax, boundary, 'r', 'old')
##fn.plot_observed(observed_new, resonance, fig, ax, boundary, 'b', 'new')
fn.plot_observed(observed_control, resonance, fig, ax, boundary, 'b', 'new')
plt.show()
