import matplotlib.pyplot as plt
import argparse

import func as fn

#resonance = '32'
#completed_path = '../../../sim_results/incl=10_completed'

parser = argparse.ArgumentParser()
parser.add_argument('--resonance', '-res',
        dest='resonance',
        action='store')
parser.add_argument('--results_path', '-results',
        dest='results',
        action='store')
args = parser.parse_args()

res_str = args.resonance
results_path = args.results
sim_results = fn.MM_sim_results(results_path, res_str)

fig, ax, boundary = fn.stability_fig_setup(res_str, results_path)
fn.plot_sims(sim_results, fig, ax)
coords = fn.interactive_mu1mu2(results_path, res_str, sim_results, fig)
plt.show()
