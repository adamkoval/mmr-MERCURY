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
from matplotlib.lines import Line2D


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


def read_biginfo(completed_path, res_str, bigin, infoout):
    """
    Reads big.in and info.out to obtain sim info.
    In:
        > completed_path - (str) path to directory
        containing completed simulations
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
        > bigin - (str) current big.in iteration, e.g.,
        0-big.in
        > infoout - (str) current info.out iteration,
        e.g., 99-info.out
    Out:
        > curr_sim - (dict) information about the
        current iteration of the simulation.
    """
    global AU, Msol
    big = '{}/{}/input/{}'.format(completed_path, res_str, bigin)
    info = '{}/{}/info/{}'.format(completed_path, res_str, infoout)
    status = get_status(info)
    with open(big) as f:
        lines = f.readlines()
        curr_sim = {'name': '{}_{}'.format(res_str, big),
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


def MM_sim_results(completed_path, res_str):
    """
    Creates list of results of all iterations of simulations
    for the resonance under consideration.
    In:
        > completed_path - (str) path to directory
        containing completed simulations
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
    Out:
        > sim_results - (1xN list) info on outcomes of
        all iterations of simulation for resonance under
        consideration.
    """
    sim_results = []
    bigins = listdir_nohidden('{}/{}/input/'.format(completed_path, res_str))
    infoouts = listdir_nohidden('{}/{}/info/'.format(completed_path, res_str))
    print(' ~~~~~~~~~~~~~~~~~~~~~~~~\n',
          'func.py/MM_sim_results():\n')
    for i in range(len(bigins)):
        sim_results.append(read_biginfo(completed_path, res_str, bigins[i], infoouts[i]))
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
    global G
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


def get_timeevol_data(completed_path, res_str, sim_results, mu1, mu2):
    """
    Get simulation data of both planets and the analytically
    determined path of planet 2.
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
        > planet1, planet2 - (dictionaries) simulation data
        of both planets. Columns are: 'Time (years)', 'a',
        'e', 'i', 'peri', 'node', 'M', 'mass'
        > model_planet2 - (dictionary) analytical migration
        data. Columns are: 'Time (years)', 'a'.
    """
    clicked_mu1 = '{:.3e}'.format(float(mu1))
    clicked_mu2 = '{:.3e}'.format(float(mu2))
    system, = [sys for sys in sim_results if '{:.3e}'.format(sys['pimass'])==clicked_mu1 and '{:.3e}'.format(sys['pomass'])==clicked_mu2]
    sim_idx = re.search('([0-9]+)-big\.in', system['name']).group(1)
    aeifile1 = "{}/{}/planets/{}-planet{}.aei".format(completed_path, res_str, sim_idx, 1)
    aeifile2 = "{}/{}/planets/{}-planet{}.aei".format(completed_path, res_str, sim_idx, 2)

    planet1 = read_planetaei(aeifile1)
    planet2 = read_planetaei(aeifile2)

    Delta = system['poDelta']
    a_fin = system['poa'] + Delta
    tau = system['potau']
    outcome = system['status'][0]
    t = planet2['Time (years)']
    model_planet2 = {'Time (years)': t, 'a': analytical_mig(t, a_fin, Delta, tau)}

    return planet1, planet2, model_planet2, sim_idx, outcome


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
    planet1, planet2, model_planet2, sim_idx, outcome = get_timeevol_data(completed_path, res_str, sim_results, mu1, mu2)
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
    ax[2].plot(t_phi, lpdiff, '.', label='$\Delta\omega$', ms=.5, alpha=_alpha)
    ax[2].axhline(0, color='k', ls='--', lw='.8', alpha=.6)
    ax[2].set_xlabel("Time [years]")
    ax[2].set_ylabel("$\Delta\phi$ [deg]")
    ax[2].legend()

    for _ax in ax:
        _ax.set_xscale('log')
        if not _ax==ax[2]:
            _ax.set_xticklabels([])
        _ax.tick_params(axis='x', which='both', direction='in')
        _ax.set_xlim(right=planet1['Time (years)'].values[-1])

    ax[0].set_title("Sim no. = {}".format(sim_idx), pad=45)
    plt.figtext(.5, .9, "$\mu_1={}$, $\mu_2={}$, res$=${}:{}, outcome={}".format(mu1, mu2, res_str[0], res_str[1], outcome), horizontalalignment='center')
    plt.tight_layout()
    plt.show()


# # # # # # # # # # # # # # #
# STABILITY BOUNDARY
# # # # # # # # # # # # # # #
def stability_boundary(res_str):
    """
    Obtain x + y data of analytical stability boundary as given by
    "Dynamics of Systems of Two Close Planets", Gladman, B., 1993.
    In:
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
    Out:
        > X, Y, mumin, mumax - see stab() below
        > x[1:], y[1:] - x and y arrays up to the x=y
        point (as stability boundary is symmetric in the
        x=y axis).
    """
    def fun(x, y, leeDelta):
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

    def stab(res_float):
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

    res_float = int(res_str[0]) / int(res_str[-1])
    X, Y, leeDelta = stab(res_float)
    y = np.linspace(0, Y[1], 100)
    x = [float(optimize.brentq(fun, 0, X[2]+.01, args=(y[i], leeDelta))) for i in range(0, len(y))]
    return X, Y, x[1:], y[1:]


# # # # # # # # # # # # # # #
# FIGURE
# # # # # # # # # # # # # # #
def stability_fig_setup(res_str):
    """
    Set up main figure displaying mu1-mu2.
    In:
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
    Out:
        > fig, ax - (objects) matplotlib.pylot figure
        and axis objects
        > boundary - (1xN array) information on stability
        boundary [see output of stability_boundary() function].
    """
    m_lims = {'21': (0, 14.15e-3), '53': (0, 4.15e-3), '32': (0, 4.5e-3), '75': (0, 2.87e-3), '43': (0, 1.42e-3)}
    fig, ax = plt.subplots()
    ax.set_xlim(m_lims[res_str])
    ax.set_ylim(m_lims[res_str])
    ax.set_xlabel('$\mu_1\ [M_1/M_\odot]$')
    ax.set_ylabel('$\mu_2\ [M_2/M_\odot]$')
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 12))
    ax.set_title('{}:{}'.format(*res_str))
    boundary = plot_boundary(res_str, fig, ax)
    return fig, ax, boundary


def plot_boundary(res_str, fig, ax):
    """
    Plots the stability boundary given by Gladman (1993) [see
    stability_boundary() function] on a given figure and axis
    object.
    In:
        > res_str - (str) specifies the resonance
        under consideration, e.g., '53', '5:3', '5-3'
        > fig, ax - (objects) matplotlib.pylot figure
        and axis objects
    Out:
        > boundary - (1xN array) information on stability
        boundary [see output of stability_boundary() function].
    """
    boundary = stability_boundary(res_str)
    xdat = np.append(boundary[2], np.flip(boundary[3]))
    ydat = np.append(boundary[3], np.flip(boundary[2]))
    ax.plot(xdat, ydat, 'b--', label='Analytical limit', linewidth=1)
    return boundary


def plot_observed(observed, res_str, fig, ax, boundary, color, label):
    """
    Plots observed systems on a given figure and axis object.
    In:
        > observed - dictionary containing all systems detected to be
        in resonance from a given online repository [see find_resonances()
        function for details]. Dictionary keys are ['21', '53', '32', '75',
        '43'], and sub-keys can be found by referring to the function
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
        > fig, ax - (objects) matplotlib.pylot figure
        and axis objects
        > boundary - (1xN array) information on stability
        boundary [see output of stability_boundary() function]
        > color - desired color of points
        > label - desired label of points
    Out:
         > fig, ax - (objects) matplotlib.pylot figure
        and axis objects.
    """
    mumin = boundary[2][-1]
    mumax = boundary[2][0]
    for system in observed[res_str]:
        pi_mass = system['pi_mass']
        sig_pi_mass = system['sig_pi_mass']
        po_mass = system['po_mass']
        sig_po_mass = system['sig_po_mass']
        s_mass = system['s_mass']
        sig_s_mass = system['sig_s_mass']
        name = system['name']
        mu_i, errs_i = mu_error(pi_mass, sig_pi_mass, s_mass, sig_s_mass)
        mu_o, errs_o = mu_error(po_mass, sig_po_mass, s_mass, sig_s_mass)
        errs_i = [[errs_i[0]], [errs_i[1]]]
        errs_o = [[errs_o[0]], [errs_o[1]]]
        ax.errorbar(mu_i, mu_o, xerr=errs_i, yerr=errs_o,
                ecolor=color, elinewidth=.5, capsize=2, fmt='.', color=color, label=label)
        first_cond = - (mumax - mumin)/mumin*mu_i + mumax
        second_cond = - mumin/(mumax - mumin)*mu_i + mumax*mumin/(mumax - mumin)
        try:
            if mu_o > first_cond and mu_o > second_cond:
                ax.annotate(name, (mu_i, mu_o))
        except:
            pass
    return fig, ax


def plot_sims(sim_results, fig, ax):
    """
    Plots simulation results on give figure and axis
    objects.
    In:
        > sim_results - (dictionary) information on results
        of simulation of the resonance under consideration
        > fig, ax - (objects) matplotlib.pylot figure
        and axis objects
    Out:
        > (No output) - fig + ax objects are affected.
    """
    print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n",
          "func.py/plot_sims():\n")
    for i, sim in enumerate(sim_results):
        x = sim['pimass'] / sim['smass']
        y = sim['pomass'] / sim['smass']
        status = sim['status'][0]
        _picker = 4.5
        if status == 'stable':
            ax.plot(x, y, 'k.', ms=8, mew=.8, fillstyle='none', picker=_picker)
        elif status == 'hit star':
            ax.plot(x, y, marker='*', c=((1, .7, .2)), ms=7, mew=.8, fillstyle='none', picker=_picker)
        elif status == 'hit planet':
            ax.plot(x, y, marker='.', c=((1, 0, 0)), ms=9, mew=.9, fillstyle='none', picker=_picker)
        elif status == 'ejected':
            ax.plot(x, y, marker='^', c=((1, .4, .75)), ms=6, mew=.8, fillstyle='none', picker=_picker)
        elif status == 'empty':
            ax.plot(x, y, marker='s', c=((0, 0, 1)), ms=5, mew=.9, fillstyle='none', picker=_picker)
        else:
            pass

        if i%100 == 0:
            print(" sim: {}".format(i))
    print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n")

    _handles = [Line2D([], [], color='k', marker='.', ms=8, mew=.8,
                ls='', fillstyle='none', label='Stable'),
            Line2D([], [], color=((1, .7, .2)), marker='*', ms=7, mew=.8,
                ls='', fillstyle='none', label='P-* collision'),
            Line2D([], [], color=((1, 0, 0)), marker='.', ms=9, mew=.9,
                ls='', fillstyle='none', label='P-P collision'),
            Line2D([], [], color=((1, .4, .75)), marker='^', ms=6, mew=.8,
                ls='', fillstyle='none', label='Ejection'),
            Line2D([], [], color=((0, 0, 1)), marker='s', ms=5, mew=.9,
                ls='', fillstyle='none', label='No data')]
    plt.legend(loc=2, handles=_handles, fancybox=True, prop={'size': 6})
    fig.subplots_adjust(bottom=.15)
    ax.text(0.01, -0.15, "N = {}".format(len(sim_results)), transform=ax.transAxes)


def interactive_mu1mu2(completed_path, res_str, sim_results, fig):
    """
    Allows click-interaction with mu1-mu2 figure to allow the inspection of
    the time-evolution of any selected system.
    In:
        > completed_path - (str) path to directory
        containing completed simulations
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
        > sim_results - (dictionary) information on results
        of simulation of the resonance under consideration
        > fig, ax - (objects) matplotlib.pylot figure
        and axis objects (mu1-mu2 graph)
    Out:
        > (No output) - time-evolution graph is plotted
        (see plot_timeevol() function).
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
        plot_timeevol(completed_path, res_str, sim_results, mu1, mu2)
        print(" ~~~~~~~~~~~~~~~~~~~~~~~~\n")
    global coords
    coords = []
    cid = fig.canvas.mpl_connect('pick_event', on_pick)
