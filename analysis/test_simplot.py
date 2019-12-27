import matplotlib.pyplot as plt

import funcs as fn

resonance = '53'
completed_path = '../completed'
sim_results = fn.MM_sim_results(completed_path, resonance)

fig, ax, boundary = fn.stability_fig_setup(resonance)
fn.plot_sims(sim_results, resonance, fig, ax)
plt.show()
