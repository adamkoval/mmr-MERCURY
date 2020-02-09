import matplotlib.pyplot as plt

import func as fn

#resonance = '32'
#completed_path = '../../../sim_results/incl=10_completed'
resonance = '21'
completed_path = '../completed'
sim_results = fn.MM_sim_results(completed_path, resonance)

fig, ax, boundary = fn.stability_fig_setup(resonance)
fn.plot_sims(sim_results, fig, ax)
coords = fn.interactive_mu1mu2(completed_path, resonance, sim_results, fig)
plt.show()
