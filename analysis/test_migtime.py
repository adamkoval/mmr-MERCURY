import numpy as np
import matplotlib.pyplot as plt

import func as fn

time = np.linspace(0, 1e7, 1000)
a_fin = 5
Delta = -3

fig, ax = plt.subplots()
for _tau in [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e15]:
    path = fn.analytical_mig(time, a_fin, Delta, _tau)
    ax.plot(time, path, label='{:.0e}'.format(_tau))

#ax.set_xscale('log')
plt.legend()
plt.show()
