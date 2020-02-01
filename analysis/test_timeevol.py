import matplotlib.pyplot as plt

import func as fn

res_str = '21'
completed_path = '../completed'
sim_results = fn.MM_sim_results(completed_path, res_str)
mu1, mu2 = 0.0116, 0.0125

fn.plot_timeevol(res_str, sim_results, mu1, mu2)
plt.show()
