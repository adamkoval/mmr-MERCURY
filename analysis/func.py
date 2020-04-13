##
#
# Functions used for finding resonances from online repository
#
##
import numpy as np
import sympy as sym
import pandas as pd
import matplotlib.pyplot as plt
import re
import os
from scipy import optimize
from scipy.signal import periodogram
from matplotlib.lines import Line2D
from matplotlib.cm import (get_cmap, ScalarMappable)
from adjustText import adjust_text

Msol = 1.9886e30 #kg
Mjup = 1.89813e27 #kg
G = 6.6741e-11 #m3kg-1s-2
AU = 1.496e11 #m


# # # # # # # # # # # # # # #
# GENERAL
# # # # # # # # # # # # # # #
def read_pathfile(pathfile):
    """
    Reads in pathfile.txt to assign paths for
    the script to use.
    In:
        > pathfile - (str) path to pathfile.txt
    Out:
        > completed - (str) path to completed
        directory.
    """
    with open(pathfile, 'r') as f:
        lines = [line for line in f.readlines() if line[0] is not '#']
        names = []
        for line in lines:
            name = line.split()[0]
            value = line.split()[2]
            globals()[name] = value
            names.append(name)
        return [globals()[name] for name in names]


def mu_error(pmass, pmasserrs, smass, smasserrs):
    """
    Function to return planet/star mass ratio from
    given planet and star masses and its error.
    In:
        > pmass - (float) planet mass in Mjup
        > pmasserrs - (1x2 tuple) errors in planet mass
        of form (min, max)
        > smass - (float) star mass in Msol
        > smasserrs - (1x2 tuple) errors in star mass
    Out:
        > mu - (float) pmass/smass mass ratio
        > (sig_minus, sig_plus) - (1x2 tuple) errors in mass
        ratio.
    """
    global Msol, Mjup
    mu = float(pmass) / float(smass) * Mjup / Msol
    sig_pmass_plus = abs((float(pmass)+float(pmasserrs[1])) / float(smass) * Mjup / Msol - mu)
    sig_pmass_minus = abs((float(pmass)+float(pmasserrs[0])) / float(smass) * Mjup / Msol - mu)
    sig_smass_plus = abs(float(pmass) / (float(smass)+float(smasserrs[1])) * Mjup / Msol - mu)
    sig_smass_minus = abs(float(pmass) / (float(smass)+float(smasserrs[0])) * Mjup / Msol - mu)
    sig_plus = np.sqrt(sig_pmass_plus**2 + sig_smass_plus**2)
    sig_minus = np.sqrt(sig_pmass_minus**2 + sig_smass_minus**2)
    return mu, (sig_minus, sig_plus)


def listdir_nohidden(path):
    """
    Function to list everything in a directory which isn't hidden,
    e.g., .swp files.
    In:
        > path - (str) path to directory to list
    Out:
        > dirlist - (1xN list) listed contents of directory.
    """
    dirlist = [thing for thing in os.listdir(path) if not thing.startswith('.')]
    return dirlist


def snip_path(path):
    """Adapts path string to format used by all functions.
    In:
        > path - (str) path to directory
    Out:
        > path - (str) amended path.
    """
    if path[-1] == "/":
        path = path[:-1]
    return path


def read_initconds(completed_path, var_str):
    """
    Picks out variable to read from initial conditions file supplied
    in each results directory.
    In:
        > completed_path - (str) path to directory
        containing completed simulations
        > var_str - (str) string of desired variable,
        as written in init_conds.txt file
    Out:
        > var - (str) variable from init_conds.txt file.
    """
    with open("{}/init_conds.txt".format(completed_path), 'r') as f:
        lines = f.readlines()
    var, = [line.split()[2] for line in lines if line.startswith(var_str)]
    return var


def get_dadt(a_f, a_t, tau):
    """Calculates radial migration speed"""
    return (a_f - a_t) / tau


# # # # # # # # # # # # # # #
# RESONANCE FINDING FUNCTIONS
# # # # # # # # # # # # # # #
def find_resonances(planets, tol):
    """
    Function for finding resonances from database provided by
    "https://exoplanetarchive.ipac.caltech.edu/".
    In:
        > planets - pandas dataframe
        > tolerance - (float) or (int) type between 0-100
                      (i.e., percentage tolerance)
        > output - (string) type specifying the desired output
                   from 'pmasses', 'smasses', 'names', 'mtypes' or 'merrors'
                   (NOTE: if 'output' is left unspecified, all values are
                   returned)
    Out:
        NOTE: all are (arrays) of (float) type values, and output depends
              on the setting of 'output' argument above.
        > pmasses - planet masses
        > smasses - star masses
        > names - system names
        > mtypes - mass types. i.e., in Mass or Msin(i) units
        > merrors - errors in mass
        > OR all of the above (see 'output' argument above)
    """
    res_floats = [2/1, 5/3, 3/2, 7/5, 4/3]
    res_strs = ['21', '53', '32', '75', '43']
    g = 0
    systems = {}
    for i in range(1, len(planets)):
        if planets['pl_hostname'][i] == planets['pl_hostname'][g] and i != len(planets)-1:
            continue
        elif planets['pl_hostname'][i] != planets['pl_hostname'][g] or i == len(planets)-1:
            sys_name = planets['pl_hostname'][g]
            if i == len(planets)-1:
                systems[sys_name] = planets[:][g:i+1]
            else:
                systems[sys_name] = planets[:][g:i]
            g = i

    obs_data = {'21': [], '53': [], '32': [], '75': [], '43': []}
    sys_names = [system for system in systems]

    for i in range(len(systems)):
        curr_sys = systems[sys_names[i]]
        for j in range(len(curr_sys)-1):
            for k in range(j+1, len(curr_sys)):
                massj = np.asarray(curr_sys['pl_bmassj'])[j]
                massk = np.asarray(curr_sys['pl_bmassj'])[k]
                if str(massj) == 'nan' or str(massk) == 'nan':
                    pass
                else:
                    periodj = np.asarray(curr_sys['pl_orbper'])[j]
                    periodk = np.asarray(curr_sys['pl_orbper'])[k]
                    if periodj > periodk:
                        outer = j
                        inner = k
                    else:
                        outer = k
                        inner = j
                    p_i_minerr = np.asarray(curr_sys['pl_orbpererr2'])[inner]
                    p_i = np.asarray(curr_sys['pl_orbper'])[inner]
                    p_i_maxerr = np.asarray(curr_sys['pl_orbpererr1'])[inner]
                    p_o_minerr = np.asarray(curr_sys['pl_orbpererr2'])[outer]
                    p_o = np.asarray(curr_sys['pl_orbper'])[outer]
                    p_o_maxerr = np.asarray(curr_sys['pl_orbpererr1'])[outer]

                    errs = {'21': [],'53': [], '32': [], '75': [], '43': []}
                    for l, res_float in enumerate(res_floats):
                        """
                        Checking the period ratios against the tolerance
                        limit for each resonance.
                        """
                        res_str = res_strs[l]

                        max_periods = [p_o+p_o_maxerr, p_i+p_i_maxerr]
                        min_periods = [p_o+p_o_minerr, p_i+p_i_minerr]
                        max_ratio = max_periods[0] / min_periods[1]
                        min_ratio = min_periods[0] / max_periods[1]
                        max_per_diff = max_ratio / res_float
                        min_per_diff = min_ratio / res_float
                        upp_per_err = abs(1 - max_per_diff)
                        low_per_err = abs(1 - min_per_diff)
                        if upp_per_err > low_per_err:
                            errs[res_strs[l]] = upp_per_err * 100
                        else:
                            errs[res_strs[l]] = low_per_err * 100
                    if min(errs.values()) < tol:
                        """
                        If the system period ratios are inside the tolerance limit,
                        append their details to the main dictionary.
                        """
                        curr_sys_info = {'name': [],
                                         'pi_mass': [],
                                         'sig_pi_mass': [],
                                         'po_mass': [],
                                         'sig_po_mass': [],
                                         's_mass': [],
                                         'sig_s_mass': [],
                                         'm_type': [],
                                         'pi_per': [],
                                         'sig_pi_per': [],
                                         'po_per': [],
                                         'sig_po_per': []}
                        curr_sys_info['name'] = '{} {}, {}'.format(
                                np.asarray(curr_sys['pl_hostname'])[inner],
                                np.asarray(curr_sys['pl_letter'])[inner],
                                np.asarray(curr_sys['pl_letter'])[outer])
                        curr_sys_info['pi_mass'] = np.asarray(curr_sys['pl_bmassj'])[inner]
                        curr_sys_info['sig_pi_mass'] = (
                                np.asarray(curr_sys['pl_bmassjerr2'])[inner],
                                np.asarray(curr_sys['pl_bmassjerr1'])[inner])
                        curr_sys_info['po_mass'] = np.asarray(curr_sys['pl_bmassj'])[outer]
                        curr_sys_info['sig_po_mass'] = (
                                np.asarray(curr_sys['pl_bmassjerr2'])[outer],
                                np.asarray(curr_sys['pl_bmassjerr1'])[outer])
                        curr_sys_info['s_mass'] = np.asarray(curr_sys['st_mass'])[inner]
                        curr_sys_info['sig_s_mass'] = (
                                np.asarray(curr_sys['st_masserr2'])[inner],
                                np.asarray(curr_sys['st_masserr1'])[inner])
                        curr_sys_info['m_type'] = np.asarray(curr_sys['pl_bmassprov'])[inner]
                        curr_sys_info['pi_per'] = np.asarray(curr_sys['pl_orbper'])[inner]
                        curr_sys_info['sig_pi_per'] = (
                                np.asarray(curr_sys['pl_orbpererr2'])[inner],
                                np.asarray(curr_sys['pl_orbpererr1'])[inner])
                        curr_sys_info['po_per'] = np.asarray(curr_sys['pl_orbper'])[outer]
                        curr_sys_info['sig_po_per'] = (
                                np.asarray(curr_sys['pl_orbpererr2'])[outer],
                                np.asarray(curr_sys['pl_orbpererr1'])[outer])
                        best_resonance = min(errs, key=errs.get)
                        obs_data[best_resonance].append(curr_sys_info)
                    else:
                        pass
    return obs_data


# # # # # # # # # # # # # # #
# SIM RESULTS
# # # # # # # # # # # # # # #
def get_status(info_file):
    """
    Get outcome of simulation by reading info.out file.
    Possible outcomes (quoted exactly as displayed in
    info.out file):
        >
           Integration complete.
        > planet1  was hit by planet2  at        3919.194 years
        > planet2  collided with the central body at     253003.2085643 years
        > planet2  ejected at     187908.8200906 years
    In:
        > info_file - (str) path to info.out file
    Out:
        > status - (str) outcome of simulation
        > planet - (str) which planet became unstable,
        if any
        > time - (float) timestamp of instability, in
        years.
    """
    def assign_vars(idx):
        """
        Returns planet and time variables read from a
        line at row idx.
        In:
            > idx - (int) specifies row number in info.out
            file (starting from 0)
        Out:
            > planet - (str) letter index of planet
            > time - (float) time of interaction.
        """
        planet = lines[idx].split()[0]
        time = float(lines[idx].split()[-2])
        return planet, time

    with open(info_file) as f:
        lines = f.readlines()
        try:
            idx = 41
            if 'collided with the central body' in lines[idx]:
                status, planet, time = ['hit star', *assign_vars(idx)]
            elif 'was hit by' in lines[idx]:
                status, planet, time = ['hit planet', *assign_vars(idx)]
            elif 'ejected at' in lines[idx]:
                status, planet, time = ['ejected', *assign_vars(idx)]
            elif 'Integration complete.' in lines[42]:
                status = 'stable'
                planet = np.nan
                time = float(lines[7].split()[-1])/365
            elif 'Continuing integration from dump files at' in lines[42]:
                idx = 45
                if 'collided with the central body' in lines[idx]:
                    status, planet, time = ['hit star', *assign_vars(idx)]
                elif 'was hit by' in lines[idx]:
                    status, planet, time = ['hit planet', *assign_vars(idx)]
                elif 'ejected at' in lines[idx]:
                    status, planet, time = ['ejected', *assign_vars(idx)]
                elif 'Integration complete.' in lines[46]:
                    status = 'stable'
                    planet = np.nan
                    time = float(lines[7].split()[-1])/365
            else:
                print(' ~~~~~~~~~~~~~~~~~~~~~~~~\n',
                      'func.py/get_status():\n',
                      'Something went wrong! Check the above function. Quitting.\n',
                      '~~~~~~~~~~~~~~~~~~~~~~~~\n')
                sys.exit()

        except IndexError:
            status = 'empty'
            planet = np.nan
            time = np.nan
    return status, planet, time


def read_biginfo(completed_path, bigin, infoout):
    """
    Reads big.in and info.out to obtain sim info.
    In:
        > completed_path - (str) path to directory
        containing completed simulations
        > bigin - (str) current big.in iteration, e.g.,
        0-big.in
        > infoout - (str) current info.out iteration,
        e.g., 99-info.out
    Out:
        > curr_sim - (dict) information about the
        current iteration of the simulation.
    """
    global AU, Msol
    completed_path = snip_path(completed_path)
    big = '{}/input/{}'.format(completed_path, bigin)
    info = '{}/info/{}'.format(completed_path, infoout)
    status = get_status(info)
    with open(big) as f:
        lines = f.readlines()
        curr_sim = {'name': '{}_{}'.format(completed_path, big),
                    'pimass': float(lines[6].split()[1][2:]), # Mjup
                    'pomass': float(lines[10].split()[1][2:]), # Mjup
                    'smass': float(1), # Msol
                    'piper': kepler3_period(float(lines[7].split()[0])*AU, Msol) / (60*60*24*365), # yrs
                    'poper': kepler3_period(float(lines[11].split()[0])*AU, Msol) / (60*60*24*365), # yrs
                    'pia': float(lines[7].split()[0]),
                    'poa': float(lines[11].split()[0]),
                    'piDelta': float(lines[6].split()[5][3:]), # AU
                    'poDelta': float(lines[10].split()[5][3:]), # AU
                    'pitau': float(lines[6].split()[4][3:]), # yrs
                    'potau': float(lines[10].split()[4][3:]), # yrs
                    'status': status
                    }
    return curr_sim


def MM_sim_results(completed_path):
    """
    Creates list of results of all iterations of simulations
    for the resonance under consideration.
    In:
        > completed_path - (str) path to directory
        containing completed simulations
    Out:
        > sim_results - (1xN list) info on outcomes of
        all iterations of simulation for resonance under
        consideration.
    """
    sim_results = []
    completed_path = snip_path(completed_path)
    bigins = listdir_nohidden('{}/input/'.format(completed_path))
    infoouts = listdir_nohidden('{}/info/'.format(completed_path))
    print(' ~~~~~~~~~~~~~~~~~~~~~~~~\n',
          'func.py/MM_sim_results():\n')
    for i in range(len(bigins)):
        sim_results.append(read_biginfo(completed_path, bigins[i], infoouts[i]))
        if i%100 == 0:
            print(' sim: {}'.format(i))
    print(' ~~~~~~~~~~~~~~~~~~~~~~~~\n')
    return sim_results


# # # # # # # # # # # # # # #
# SIM TIME EVOL
# # # # # # # # # # # # # # #
def read_planetaei(aeifile):
    """
    Reads .aei file containing time-evolution of planets.
    In:
        > aeifile - (str) path to planet.aei file
    Out:
        > planet - (pandas dict) read time-dependent
        properties of planet in simulation.
    """
    headers = ['Time (years)', 'a', 'e', 'i', 'peri', 'node', 'M', 'mass']
    planet = pd.read_csv(aeifile, skiprows=4, names=headers, delim_whitespace=True, dtype=np.float64, na_values="********")
    return planet


def kepler3_period(a, M_star):
    """
    Kepler's third law for finding orbital period from
    semi-major axis and star mass.
    In:
        > a - (float or int) semi-major axis
        > M_star - (float or int) star mass
    Out:
        > T - (float) orbital period of planetary
        companion.
    """
    global G
    T = 2*np.pi*np.sqrt(a**3/(G*M_star))
    return T


def kepler3_resdisp(res_float, a_i):
    """
    Kepler's third law for finding semi-major
    axis of outer planet from nominal resonance
    and inner semi-major axis.
    In:
        > res_float - (float) fractional value of resonance under
        consideration, i.e., 2/1=2.0, 5/3=1.66667, and so on
        > a_i - (float) semi-major axis of inner planet [units=AU]
    Out:
        > a_o - (float) semi-major axis of outer planet [units=AU].
    """
    a_o = res_float**(2/3) * a_i
    return a_o


def get_resvar(res_str, planet_i, planet_o):
    """
    In:
        > planet_i, planet_o - (csv) data read from
        planet1.aei and planet2.aei files for the
        the inner and outer planets resp.
    Out:
        > phi1 - (1xN array) time-evolution of
        resonance variable with
        {j1, j2, j3, j4} = {p, -q, -(p-1), 0}
        > phi2 - (1xN array) time-evolution of
        resonance variable with
        {j1, j2, j3, j4} = {p, -q, 0, -(p-1)}
    """
    def truncate_longer(arr1, arr2):
        """
        Truncate the longer array to the length of the
        shorter one.
        In:
            > arr1, arr2 - (1xN arrays)
        Out:
            > arr1, arr2 - (1xN arrays) after being passed
            through condition.
        """
        if len(arr1) > len(arr2):
            print(" Truncating arr1.")
            return arr1[:len(arr2)], arr2
        elif len(arr2) > len(arr1):
            print(" Truncating arr2.")
            return arr1, arr2[:len(arr1)]
        else:
            print(" Not truncating, arrays are of same length.")
            return arr1, arr2

    print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n",
          "func.py/get_resvar():\n")
    planet_i, planet_o = truncate_longer(planet_i, planet_o)
    print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n")

    omega_i = planet_i['node'] + planet_i['peri'] # Longitude of periapsis
    omega_o = planet_o['node'] + planet_o['peri']
    Lambda_i = omega_i + planet_i['M'] # Mean longitude
    Lambda_o = omega_o + planet_o['M']

    p = int(res_str[0])
    q = int(res_str[1])
    phi1 = p*Lambda_o - q*Lambda_i - (p - q)*omega_o
    phi2 = p*Lambda_o - q*Lambda_i - (p - q)*omega_i
    deltaphi = phi1 - phi2
    t_phi = planet_i['Time (years)']
    lpdiff = omega_o - omega_i
    return phi1, phi2, t_phi, deltaphi, lpdiff


def get_timeevol_data(completed_path, sim_results, mu1, mu2):
    """
    Get simulation data of both planets and the analytically
    determined path of planet 2.
    In:
        > completed_path - (str) path to directory
        containing completed simulations
        > sim_results - (dictionary) information on results
        of simulation of the resonance under consideration
        > mu1, mu2 - (float) mass ratios of inner and outer
        planet (resp.), to be passed to get_timeevol_data()
        fucntion
    Out:
        > planet1, planet2 - (dictionaries) simulation data
        of both planets. Columns are: 'Time (years)', 'a',
        'e', 'i', 'peri', 'node', 'M', 'mass'
        > model_planet2 - (dictionary) analytical migration
        data. Columns are: 'Time (years)', 'a'.
    """
    completed_path = snip_path(completed_path)
    clicked_mu1 = '{:.3e}'.format(float(mu1))
    clicked_mu2 = '{:.3e}'.format(float(mu2))

    system, = [sys for sys in sim_results if '{:.3e}'.format(sys['pimass'])==clicked_mu1 and '{:.3e}'.format(sys['pomass'])==clicked_mu2]
    sim_idx = re.search('([0-9]+)-big\.in', system['name']).group(1)
    aeifile1 = "{}/planets/{}-planet{}.aei".format(completed_path, sim_idx, 1)
    aeifile2 = "{}/planets/{}-planet{}.aei".format(completed_path, sim_idx, 2)

    planet1 = read_planetaei(aeifile1)
    planet2 = read_planetaei(aeifile2)

    Delta = system['poDelta']
    a_fin = system['poa'] + Delta
    tau = system['potau']
    outcome, planet, final_age = system['status']
    t = planet2['Time (years)']

    model_planet2 = {'Time (years)': t, 'a': analytical_mig(t, a_fin, Delta, tau)}

    return planet1, planet2, model_planet2, sim_idx, tau, outcome, final_age, planet


def analytical_mig(t, a_fin, Delta, tau):
    """Analytical migration model by Malhotra, R., ``The
    origin of Pluto's peculiar orbit", Letters to Nature, 1993.
    In:
        > t - (1xN array) time array from planet.aei
        file
        > a_fin - (float) final semi-major axis after
        migration [unit=AU]
        > Delta - (float) migration distance [unit=AU]
        > tau - (float) migration timescale [unit=s]
    Out:
        > a_t - (1xN array) semi-major axis evolution
        of planet2 [unit=AU].
    """
    a_t = [(a_fin - Delta*np.exp(-_t/tau)) for _t in t]
    return a_t


def plot_timeevol(completed_path, res_str, sim_results, mu1, mu2):
    """
    Plot time-evolution graph from simulation data.
    (NOTE: makes use of get_timeevol_data() function)
    In:
        > completed_path - (str) path to directory
        containing completed simulations
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
        > sim_results - (dictionary) information on results
        of simulation of the resonance under consideration
        > mu1, mu2 - (float) mass ratios of inner and outer
        planet (resp.), to be passed to get_timeevol_data()
        fucntion
    Out:
        > (No output) time-evolution graph is plotted and
        displayed.
    """
    completed_path = snip_path(completed_path)
    planet1, planet2, model_planet2, sim_idx, tau, outcome, final_age, planet = get_timeevol_data(completed_path, sim_results, mu1, mu2)
    phi1, phi2, t_phi, deltaphi, lpdiff = get_resvar(res_str, planet1, planet2)
    res_float = float(res_str[0])/float(res_str[1])
    a_i = planet1['a'][0]
    a_o = kepler3_resdisp(res_float, a_i)
#    system, = [sys for sys in sim_results if sim_idx==re.search('([0-9]+)-big\.in', system['name']).group(1)]

    fig, ax = plt.subplots(3)
    linewidth = .6
    _alpha = .76

    ax[0].plot(planet1['Time (years)'], planet1['a'], label='planet 1', lw=linewidth, alpha=_alpha)
    ax[0].plot(planet2['Time (years)'], planet2['a'], label='planet 2', lw=linewidth, alpha=_alpha)
    ax[0].plot(model_planet2['Time (years)'], model_planet2['a'], label='planet 2 (model)')
    ax[0].axhspan(a_i, a_o, alpha=.14, color='k', label='{}:{} resonance'.format(*res_str))
    ax[0].set_ylabel("a [A. U.]")
    ax[0].legend(bbox_to_anchor=(0, .985, 1, 0), loc=3, ncol=4, mode='expand',
            fancybox=True, prop={'size': 9})

    ax[1].plot(planet1['Time (years)'], planet1['e'], label='planet1', lw=linewidth, alpha=_alpha)
    ax[1].plot(planet2['Time (years)'], planet2['e'], label='planet2', lw=linewidth, alpha=_alpha)
    ax[1].set_ylabel("e [no unit]")

    #ax[2].plot(t_phi, phi1, label='phi1', lw=linewidth, alpha=_alpha)
    #ax[2].plot(t_phi, phi2, label='phi2', lw=linewidth, alpha=_alpha)
    #ax[2].plot(t_phi, deltaphi, label='Resonance variable difference', c='r', lw=linewidth, alpha=_alpha)
    ax[2].plot(t_phi, lpdiff, '.', ms=.5, alpha=_alpha)
    ax[2].axhline(0, color='k', ls='--', lw='.8', alpha=.6)
    ax[2].set_xlabel("Time [years]")
    ax[2].set_ylabel("$\\Delta\\varpi$ [deg]")

    for _ax in ax:
        _ax.set_xscale('log')
        if not _ax==ax[2]:
            _ax.tick_params(labelbottom=False)
        _ax.tick_params(axis='x', which='both', direction='in')
        #_ax.set_xlim(right=planet1['Time (years)'].values[-1])
        _ax.set_xlim(right=1e7)


    ax[0].set_title("Sim no. = {}".format(sim_idx), pad=45)
    plt.figtext(.5, .9,
            "$\mu_1={}$, $\mu_2={}$, res$=${}:{}, outcome={}, planet={}, tau={:.2e} (years)".format(mu1, mu2, res_str[0], res_str[1], outcome, planet, tau),
            horizontalalignment='center')
    pratio = get_pratio(planet1['a'], planet2['a'])
    plt.figtext(.5, .01,
            "Final period ratio: {:.4g}".format(pratio),
            )
    plt.tight_layout()
    plt.show()


def plot_periods(completed_path, res_str, sim_results, mu1, mu2):
    """
    Plot period analysis graph from simulation data.
    (NOTE: makes use of get_timeevol_data() function)
    In:
        > completed_path - (str) path to directory
        containing completed simulations
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
        > sim_results - (dictionary) information on results
        of simulation of the resonance under consideration
        > mu1, mu2 - (float) mass ratios of inner and outer
        planet (resp.), to be passed to get_timeevol_data()
        fucntion
    Out:
        > (No output) time-evolution graph is plotted and
        displayed.
    """
    completed_path = snip_path(completed_path)
    planet1, planet2, model_planet2, sim_idx, tau, outcome, final_age, planet = get_timeevol_data(completed_path, sim_results, mu1, mu2)
    phi1, phi2, t_phi, deltaphi, lpdiff = get_resvar(res_str, planet1, planet2)

    fig, ax = plt.subplots(4)
    linewidth = .6
    _alpha = .76

    # For period analysis
    fs = 1/(1000*365*24*60*60)
    t = planet1['Time (years)'][1000:]
    a_p1 = planet1['a'][1000:]
    a_p2 = planet2['a'][1000:]
    e_p1 = planet1['e'][1000:]
    e_p2 = planet2['e'][1000:]
    lpd = lpdiff[1000:]

    a_p1_f, a_p1_Pxx = periodogram(a_p1, fs)
    a_p2_f, a_p2_Pxx = periodogram(a_p2, fs)
    ax[0].plot(a_p1_f, a_p1_Pxx, label="p1 a frequencies")
    ax[1].plot(a_p2_f, a_p2_Pxx, label="p2 a frequencies")

    e_p1_f, e_p1_Pxx = periodogram(e_p1, fs)
    e_p2_f, e_p2_Pxx = periodogram(e_p2, fs)
    ax[2].plot(e_p1_f, e_p1_Pxx, label="p1 e frequencies")
    ax[3].plot(e_p2_f, e_p2_Pxx, label="p2 e frequencies")

    for _ax in ax:
        if not _ax==ax[3]:
            _ax.tick_params(labelbottom=False)
        _ax.tick_params(axis='x', which='both', direction='in')
        _ax.legend()

    plt.tight_layout()


def get_pratio(a1, a2):
    """Finds final period ratio of stable systems by averaging over last
    1e6 years.
    """
    global Msol
    final_avg_a1 = np.mean(a1[1000:])
    final_avg_a2 = np.mean(a2[1000:])
    p1 = kepler3_period(final_avg_a1, Msol)
    p2 = kepler3_period(final_avg_a2, Msol)
    pratio = p2/p1
    return pratio

# # # # # # # # # # # # # # #
# STABILITY BOUNDARY
# # # # # # # # # # # # # # #
class StabilityBoundary:
    def __init__(self, res_str):
        self.res_str = res_str
        self.res_float = int(res_str[0]) / int(res_str[-1])
        X, Y, leeDelta = self.stab(self.res_float)
        y = np.linspace(0, Y[1], 100)
        x = [float(optimize.brentq(self.fun, 0, X[2]+.01, args=(y[i], leeDelta))) for i in range(0, len(y))]
        self.x = x[1:]
        self.y = y[1:]


    def fun(self, x, y, leeDelta):
        """
        Function for fitting, as given by Gladman (1993). To be
        used for finding stability boundary.
        In:
            > x, y - (undefined) independent variables
            > leeDelta - (float) normalised separation of semi-major
            axes of planet pair
        Out:
            > eqn23 - (function) of the form of eqn. 23 in
            the source paper.
        """
        eqn23 = 2*3**(1/6)*(x + y)**(1/3) + 2*3**(1/3)*(x + y)**(2/3) - (11*x + 7*y)/(3**(11/6)*(x + y)**(1/3)) - leeDelta
        return eqn23


    def stab(self, res_float):
        """
        Finds limits of stability boundary (i.e., minimum and
        maximum mu values up to the x=y point).
        In:
            > res_float - (float) fractional value of resonance under
        consideration, i.e., 2/1=2.0, 5/3=1.66667, and so on
        Out:
            > X, Y - (1x3 arrays) denote the maximum point,
            x=y symmetry point and minimum point, in both axes.
            NOTE: X[1]==Y[1]
            > leeDelta - (float) normalised separation of semi-major
            axes of planet pair.
        """
        x = sym.Symbol('x')
        leeDelta = res_float**(2/3) - 1
        mumax = sym.solve(2*3**(1/6) * (x)**(1/3) + 2*3**(1/3) * (x)**(2/3) - 11*x / (3**(11/6) * x**(1/3)) - leeDelta, x)
        mumin = sym.solve(2*3**(1/6) * (2*x)**(1/3) + 2*3**(1/3) * (2*x)**(2/3) - 18*x / (3**(11/6) * (2*x)**(1/3)) - leeDelta, x)
        X = [0, float(mumin[0]), float(mumax[0])]
        Y = [float(mumax[0]), float(mumin[0]), 0]
        return X, Y, leeDelta


# # # # # # # # # # # # # # #
# FIGURE
# # # # # # # # # # # # # # #
class StabilityFigure:
    def __init__(self, res_str):
        self.res_str = res_str
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlabel("$\mu_1\ [M_1/M_\odot]$")
        self.ax.set_ylabel("$\mu_2\ [M_2/M_\odot]$")
        plt.ticklabel_format(axis='both', style='sci', scilimits=(-1, 1))
        #self.ax.set_title('{}:{} ({})'.format(*res_str, completed_path)) # DEBUG
        self.ax.set_title("{}:{}".format(*res_str))
        self.boundary = StabilityBoundary(res_str)
        self.plot_boundary()
        self.handles = []

    def show(self):
        """Display plot"""
        #self.m_lims = {'21': (0, 14.15e-3), # ORIGINAL
        #        '53': (0, 4.15e-3),
        #        '32': (0, 4.5e-3),
        #        '75': (0, 2.87e-3),
        #        '43': (0, 1.42e-3)}
        #self.m_lims = {'21': (0, 20e-3), # EXPANDED
        #        '53': (0, 4.15e-3),
        #        '32': (0, 3.5e-3),
        #        '75': (0, 2.87e-3),
        #        '43': (0, 1.42e-3)}
        #lim = ylim = self.m_lims[self.res_str]
        xmax = self.ax.get_xlim()[1]
        ymax = self.ax.get_ylim()[1]
        mx = max(xmax, ymax)
        self.ax.set_xlim(0, mx)
        self.ax.set_ylim(0, mx)

        plt.legend(loc=2, handles=self.handles, fancybox=True, prop={'size': 6})
        try:
            adjust_text(self.observed_labels)
        except AttributeError:
            pass
        plt.show()


    def plot_boundary(self):
        """
        Plots the stability boundary given by Gladman (1993) [see
        stability_boundary() function] on a given figure and axis
        object.
        Out:
            > boundary - (1xN array) information on stability
            boundary [see output of stability_boundary() function].
        """
        bound_x = self.boundary.x
        bound_y = self.boundary.y
        xdat = np.append(bound_x, np.flip(bound_y))
        ydat = np.append(bound_y, np.flip(bound_x))
        self.ax.plot(xdat, ydat, 'b--', label='Analytical limit', linewidth=1)


    def plot_observed(self, observed):
        """
        Plots observed systems on a given figure and axis object.
        In:
            > observed - dictionary containing all systems detected to be
            in resonance from a given online repository [see find_resonances()
            function for details]. Dictionary keys are ['21', '53', '32', '75',
            '43'], and sub-keys can be found by referring to the function
        """
        mumin = self.boundary.x[-1]
        mumax = self.boundary.x[0]
        self.observed_labels = []
        for system in observed[self.res_str]:
            name = system['name']
            s_mass = system['s_mass']
            pi_mass = system['pi_mass']
            po_mass = system['po_mass']
            m_prov = system['m_type']
            sig_s_mass = system['sig_s_mass']
            sig_pi_mass = system['sig_pi_mass']
            sig_po_mass = system['sig_po_mass']
            mu_i, errs_i = mu_error(pi_mass, sig_pi_mass, s_mass, sig_s_mass)
            mu_o, errs_o = mu_error(po_mass, sig_po_mass, s_mass, sig_s_mass)
            errs_i = [[errs_i[0]], [errs_i[1]]]
            errs_o = [[errs_o[0]], [errs_o[1]]]

            c_Mass = (.7, .15, 1)
            c_Msini = (.1, .7, 1)
            if m_prov == "Mass" or m_prov == "Msin(i)/sin(i)":
                color = c_Mass
            elif m_prov == "Msini":
                color = c_Msini

            self.ax.errorbar(mu_i, mu_o, xerr=errs_i, yerr=errs_o,
                    ecolor=color, elinewidth=.5, capsize=2, fmt='.', ms=10, mec='k', color=color)
            first_cond = - (mumax - mumin)/mumin*mu_i + mumax
            second_cond = - mumin/(mumax - mumin)*mu_i + mumax*mumin/(mumax - mumin)
            try:
                if mu_o > first_cond and mu_o > second_cond:
                    self.observed_labels.append(self.ax.text(mu_i, mu_o, name,
                        fontsize=12, color='k', weight='bold'))
            except:
                pass
        self.handles.append(Line2D([], [], color=c_Mass, marker='.',
            ms=8, mec='k', ls='', label="Observed (Mass)"))
        self.handles.append(Line2D([], [], color=c_Msini, marker='.',
            ms=8, mec='k', ls='', label="Observed (Msini)"))


    def plot_sims(self, completed_path, view):
        """
        Plots simulation results on given figure and axis
        objects.
        In:
            > sim_results - (dictionary) information on results
            of simulation of the resonance under consideration
        """
        print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n",
              "func.py/plot_sims():\n")

        sim_results = MM_sim_results(completed_path)
        _picker = 4.5
        stable_age = 10006850
        cmap_age = get_cmap('brg')
        cmap_prat = get_cmap('jet')
        prat_ticks = [1.333, 1.4, 1.5, 1.667, 2.0]
        prat_tick_range = prat_ticks[-1] - prat_ticks[0]
        settings = {'stable': [(0, 0, 0), '.', 8, .8, "Stable"],
                'hit star': [(1, .7, .2), '*', 7, .8, "P-* collision"],
                'hit planet': [(1, 0, 0), '.', 9, .9, "P-P collision"],
                'ejected': [(1, .4, .75), '^', 6, .8, "Ejection"],
                'empty': [(0, 0, 1), 's', 5, .9, "No data"]}

        for i, sim in enumerate(sim_results):
            try:
                x = sim['pimass'] / sim['smass']
                y = sim['pomass'] / sim['smass']
                outcome, planet, final_age = sim['status']
                if view=="age":
                    fractional_age = final_age / stable_age
                    _color, _marker, _ms, _mew = cmap_age(fractional_age)[:3], '.', 8, .8
                elif view=="outcome":
                    _color, _marker, _ms, _mew = settings[outcome][:-1]
                elif view=="pratio":
                    if outcome=='stable':
                        sim_idx = re.search('([0-9]+)-big\.in', sim['name']).group(1)
                        aeifile1 = "{}/planets/{}-planet{}.aei".format(completed_path, sim_idx, 1)
                        aeifile2 = "{}/planets/{}-planet{}.aei".format(completed_path, sim_idx, 2)

                        planet1 = read_planetaei(aeifile1)
                        planet2 = read_planetaei(aeifile2)
                        a1 = planet1['a']
                        a2 = planet2['a']
                        pratio = get_pratio(a1, a2)
                        fractional_pratio = (pratio - prat_ticks[0])/prat_tick_range
                        _color = cmap_prat(fractional_pratio)[:3]
                        _marker, _ms, _mew = '.', 8, .8
                    else:
                        _color, _marker, _ms, _mew = 'k', '.', 2, .2

                self.ax.plot(x, y, color=_color, marker=_marker, ms=_ms, mew=_mew,
                        fillstyle='none', picker=_picker)

                if i%100 == 0:
                    print(" sim: {}".format(i))

            except FileNotFoundError:
                pass
        print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n")

        if view=="age":
            cbar = plt.colorbar(ScalarMappable(cmap=cmap_age), ax=self.ax)
            cbar.set_ticks([0, .2, .4, .6, .8, 1])
            cbar.set_ticklabels(['{:.1e}'.format(val)
                for val in np.dot([0, .2, .4, .6, .8, 1], stable_age)])
        elif view=="outcome":
            for outcome in settings:
                _color, _marker, _ms, _mew, _label = settings[outcome]
                _handle = Line2D([], [], color=_color, marker=_marker,
                        ms=_ms, mew=_mew, ls='', fillstyle='none',
                        label=_label)
                self.handles.append(_handle)
        elif view=="pratio":
            cbar = plt.colorbar(ScalarMappable(cmap=cmap_prat), ax=self.ax)
            cbar.set_ticks([(val-prat_ticks[0])/prat_tick_range for val in prat_ticks])
            cbar.set_ticklabels(prat_ticks)

        self.fig.subplots_adjust(bottom=.17)

        potau = sim['potau']
        poDelta = sim['poDelta']
        pia = sim['pia']
        poa = sim['poa']
        poa_f = poDelta + poa
        res_float = float(self.res_str[0])/float(self.res_str[1])
        a_res = kepler3_resdisp(res_float, pia)
        dadt = get_dadt(poa_f, a_res, potau)
        #a_res74 = kepler3_resdisp(7/4, pia) # DEBUG
        #dadt = get_dadt(poa_f, a_res74, potau) # DEBUG
        incl_2 = float(read_initconds(completed_path, "incl_2"))

        self.ax.text(0, -.2,
                "N = {}".format(len(sim_results))
                + ",\t"
                + "$\\left.\\frac{{da}}{{dt}}\\right|_{{t_{{res}}}}={:.3e}$ [AU/yr]".format(dadt)
                + ",\t"
                + "$i$ = {:.1f} [Deg]".format(incl_2),
                color=(.1, .3, .6),
                transform=self.ax.transAxes)

        self.interactive_mu1mu2(completed_path, sim_results)


    def interactive_mu1mu2(self, completed_path, sim_results):
        """
        Allows click-interaction with mu1-mu2 figure to allow the inspection of
        the time-evolution of any selected system.
        In:
            > completed_path - (str) path to directory
            containing completed simulations
            > sim_results - (dictionary) information on results
            of simulation of the resonance under consideration
        """
        def on_pick(event):
            point = event.artist
            x, y = point.get_data()
            global coords
            coords.append((x[0], y[0]))
            print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n",
                  "func/plot_sims_timeevol():\n")
            if len(coords) > 1:
                idx = len(coords) - 1
                print("You have selected more than one system.\n",
                        "If you did this by accident, please\n",
                        "zoom in to select just one at a time.\n",
                        "Plotting system {}.\n".format(coords[idx]))
                #fig.canvas.mpl_disconnect(cid)
            else:
                print("Plotting system {}.".format(coords[0]))
                idx = 0
            mu1 = coords[idx][0]
            mu2 = coords[idx][1]
            plot_timeevol(completed_path, self.res_str, sim_results, mu1, mu2)
            #plot_periods(completed_path, self.res_str, sim_results, mu1, mu2) # DEBUG
            print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n")

        #completed_path = snip_path(completed_path)
        global coords
        coords = []
        cid = self.fig.canvas.mpl_connect('pick_event', on_pick)
