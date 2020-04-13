import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

import func as fn
from matplotlib.lines import Line2D

parser = argparse.ArgumentParser()
parser.add_argument('--tolerance', '-tol',
        dest='tolerance',
        action='store',
        help="\n".join(["In units of percent. Tolerance from the nominal resonance.",
        "If absent, observed planets are not plotted."]))
parser.add_argument('--resonance', '-res_str',
        dest='res_str',
        action='store',
        help="\n".join(["Nominal resonance. Must be passed if sim results aren't",
        "plotted; otherwise redundant."]))
parser.add_argument('--results_path', '-results',
        dest='results',
        action='store',
        help="\n".join(["Relative path from this script to the results directory.",
        "e.g., ../results/completed_20-03-30.",
        "If absent, sims are not plotted."]))
parser.add_argument('--view', '-v',
        dest='view',
        action='store',
        help="\n".join(["Which format to view the mu1-mu2 graph in.",
        "Possible formats are 'age', 'outcome' or 'pratio'.",
        "Must be included if --results_path argument is passed."]))
args = parser.parse_args()

results = args.results
view = args.view
res_str = args.res_str
tol = args.tolerance

if tol and res_str and results==None:
    tol = float(tol)
    stab_fig = fn.StabilityFigure(res_str)
    planets_csv_new = pd.read_csv('planets_new.csv',
            skiprows=34, header=0, delim_whitespace=False)
    observed_new = fn.find_resonances(planets_csv_new, tol)
    stab_fig.plot_observed(observed=observed_new)
    stab_fig.show()

elif tol and results and view:
    tol = float(tol)
    res_str = fn.read_initconds(results, "res")
    res_str = res_str[0] + res_str[2]
    stab_fig = fn.StabilityFigure(res_str)
    stab_fig.plot_sims(completed_path=results, view=view)
    planets_csv_new = pd.read_csv('planets_new.csv',
            skiprows=34, header=0, delim_whitespace=False)
    observed_new = fn.find_resonances(planets_csv_new, tol)
    stab_fig.plot_observed(observed=observed_new)
    stab_fig.show()

elif results and view and tol==None:
    res_str = fn.read_initconds(results, "res")
    res_str = res_str[0] + res_str[2]
    stab_fig = fn.StabilityFigure(res_str)
    stab_fig.plot_sims(completed_path=results, view=view)
    stab_fig.show()
