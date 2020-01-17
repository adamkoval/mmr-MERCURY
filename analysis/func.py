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

Msol = 1.9886e30 #kg
Mjup = 1.89813e27 #kg
G = 6.6741e-11 #m3kg-1s-2
AU = 1.496e11 #m


# # # # # # # # # # # # # # #
# GENERAL
# # # # # # # # # # # # # # #
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
    with open(info_file) as f:
        lines = f.readlines()
        try:
            if 'collided with the central body' in lines[41]:
                status = 'hit star'
                planet = lines[41].split()[0]
                time = float(lines[41].split()[-2])
            elif 'was hit by' in lines[41]:
                status = 'hit planet'
                planet = lines[41].split()[0]
                time = float(lines[41].split()[-2])
            elif 'ejected at' in lines[41]:
                status = 'ejected'
                planet = lines[41].split()[0]
                time = float(lines[41].split()[-2])
            elif 'Integration complete.' in lines[42]:
                status = 'stable'
                planet = np.nan
                time = float(lines[7].split()[-1])/365
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
                    'piper': kepler3(float(lines[7].split()[0])*AU, Msol) / (60*60*24*365), # yrs
                    'poper': kepler3(float(lines[11].split()[0])*AU, Msol) / (60*60*24*365), # yrs
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
    for i in range(len(bigins)):
        sim_results.append(read_biginfo(completed_path, res_str, bigins[i], infoouts[i]))
    return sim_results


# # # # # # # # # # # # # # #
# SIM TIME EVOL
# # # # # # # # # # # # # # #
def read_planet(aeifile):
    """
    Reads .aei file containing time-evolution of planets.
    In:
        > aeifile - (str) path to planet.aei file
    Out:
        > planet - (pandas dict) read time-dependent
        properties of planet in simulation.
    """
    headers = ['time', 'a', 'e', 'i', 'peri', 'node', 'M', 'mass']
    planet = pd.read_csv(aeifile, skiprows=4, names=headers, delim_whitespace=True)
    return planet


def kepler3(a, M):
    """
    Kepler's third law.
    In:
        > a - (float or int) semi-major axis
        > M - (float or int) star mass
    Out:
        > T - (float) orbital period of planetary
        companion.
    """
    global G
    T = 2*np.pi*np.sqrt(a**3/(G*M))
    return T


def resvar(planet_i, planet_o, res_str):
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
            print('Truncating arr1')
            return arr1[:len(arr2)], arr2
        elif len(arr2) > len(arr1):
            print('Truncating arr2')
            return arr1, arr2[:len(arr1)]
        else:
            print('Not truncating')
            return arr1, arr2

    planet_i, planet_o = truncate_longer(planet_i, planet_o)
    Lambda_i = planet_i['node'] + planet_i['peri'] # Longidute of periapsis
    Lambda_o = planet_o['node'] + planet_o['peri']
    omega_i = Lambda_i + planet_i['M'] # Mean longitude
    omega_o = Lambda_o + planet_o['M']
    p = int(res_str[0])
    q = int(res_str[1])
    phi1 = p*Lambda_o - q*Lambda_i - (p - q)*omega_o
    phi2 = p*Lambda_o - q*Lambda_i - (p - q)*omega_i
    t_phi = planet_i['time']
    return phi1, phi2, t_phi

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


def plot_sims(sim_results, res_str, fig, ax):
    """
    Plots simulation results on give figure and axis
    objects.
    In:
        > sim_results - (dictionary) information on results
        of simulation of the resonance under consideration
        > res_str - (str) species the resonance
        under consideration, e.g., '53', '5:3', '5-3'
        > fig, ax - (objects) matplotlib.pylot figure
        and axis objects
    Out:
        > fig, ax - (objects) matplotlib.pylot figure
        and axis objects.
    """
    for sim in sim_results:
        x = sim['pimass'] / sim['smass']
        y = sim['pomass'] / sim['smass']
        status = sim['status'][0]
        if status == 'stable':
            ax.plot(x, y, 'k.', ms=3.2, mew=3.2)
        elif status == 'hit star':
            ax.plot(x, y, marker='*', c=((1, .9, .2)), ms=7, mew=.05)
        elif status == 'hit planet':
            ax.plot(x, y, marker='.', c=((1, 0, 0)), ms=1, mew=5)
        elif status == 'ejected':
            ax.plot(x, y, marker='^', c=((1, .5, 0)), ms=5, mew=.05)
        elif status == 'empty':
            ax.plot(x, y, marker='s', c=((0, 0, 1)), ms=5, mew=.05)
        else:
            pass
    return fig, ax
