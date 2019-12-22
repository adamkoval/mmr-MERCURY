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
# RESONANCE FINDING FUNCTIONS
# # # # # # # # # # # # # # #
def find_resonances(planets, tol):
    """
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
# GENERAL
# # # # # # # # # # # # # # #
def mu_error(pmass, pmasserrs, smass, smasserrs):
    global Msol, Mjup
    mu = float(pmass) / float(smass) * Mjup / Msol
    sig_pmass_plus = abs((float(pmass)+float(pmasserrs[1])) / float(smass) * Mjup / Msol - mu)
    sig_pmass_minus = abs((float(pmass)+float(pmasserrs[0])) / float(smass) * Mjup / Msol - mu)
    sig_smass_plus = abs(float(pmass) / (float(smass)+float(smasserrs[1])) * Mjup / Msol - mu)
    sig_smass_minus = abs(float(pmass) / (float(smass)+float(smasserrs[0])) * Mjup / Msol - mu)
    sig_plus = np.sqrt(sig_pmass_plus**2 + sig_smass_plus**2)
    sig_minus = np.sqrt(sig_pmass_minus**2 + sig_smass_minus**2)
    return mu, (sig_minus, sig_plus)


def mass_ratio(pmass, smass):
    global Msol, Mjup
    mu = float(pmass) / float(smass) * Mjup / Msol
    return mu


# # # # # # # # # # # # # # #
# READ IN DATA
# # # # # # # # # # # # # # #
def get_status(info_file):
    """
    Possible statuses:
        >
           Integration complete.
        > planet1  was hit by planet2  at        3919.194 years
        > planet2  collided with the central body at     253003.2085643 years
        > planet2  ejected at     187908.8200906 years
    """
    with open(info_file) as f:
        lines = f.readlines()
    if 'Integration complete.' in lines[42]:
        status = 'stable'
        planet = np.nan
        time = lines[7].split()[-1]/365
    elif 'collided with the central body' in lines[41]:
        status = 'hit star'
        planet = lines[41].split()[0]
        time = lines[41].split()[-1]
    elif 'was hit by' in lines[41]:
        status = 'hit planet'
        planet = lines[41].split()[0]
        time = lines[41].split()[-1]
    elif 'ejected at' in lines[41]:
        status = 'ejected'
        planet = lines[41].split()[0]
        time = lines[41].split()[-1]
    return status, planet, time


def read_planet(datafile):
    headers = ['time', 'a', 'e', 'i', 'peri', 'node', 'M', 'mass']
    planet = pd.read_csv(datafile, skiprows=4, names=headers, delim_whitespace=True)
    return planet

def kepler3(a, M):
    global G
    return 2*np.pi*np.sqrt(a**3/(G*M))


def read_big(path, res, run, sim):
    global AU, Msol
    datafile = '{}/{}/{}/input/{}'.format(path, res, run, sim)
    with open(datafile) as f:
        lines = f.readlines()
    curr_sim = {'name': '{}_{}_{}'.format(res, run, sim),
                'pimass': float(lines[6].split()[1][2:]), # Mjup
                'pomass': float(lines[10].split()[1][2:]),
                'smass': float(1), # Msol
                'piper': kepler3(float(lines[7].split()[0])*AU, Msol) / (60*60*24*365), # years
                'poper': kepler3(float(lines[11].split()[0])*AU, Msol) / (60*60*24*365)}
    return curr_sim


def MM_sim_data(path):
    resonances = ('2-1', '5-3', '3-2', '7-5', '4-3')
    sim_data = {'21': [], '53': [], '32': [], '75': [], '43': []}
    for res in resonances:
        try:
            runs = [run for run in os.listdir('{}/{}'.format(path, res)) if re.match('[0-9]+$', run)]
            for run in runs:
                inputs = os.listdir('{}/{}/{}/input/'.format(path, res, run))
                for bigin in inputs:
                    sim_data[res[:1] + res[-1]].append(read_big(path, res, run, bigin))
        except:
            pass
    return sim_data


# # # # # # # # # # # # # # #
# STABILITY BOUNDARY
# # # # # # # # # # # # # # #
def stability_boundary(resonance):
    def fun(x, y, delta):
        return 2 * 3**(1/6) * (x+y)**(1/3) + 2 * 3**(1/3) * (x+y)**(2/3) - (11*x+7*y) / (3**(11/6) * (x+y)**(1/3)) - delta

    def stab(resonance):
        x = sym.Symbol('x')
        leeDelta = resonance**(2/3) - 1
        mumax = sym.solve(2*3**(1/6) * (x)**(1/3) + 2*3**(1/3) * (x)**(2/3) - 11*x / (3**(11/6) * x**(1/3)) - leeDelta, x)
        mumin = sym.solve(2*3**(1/6) * (2*x)**(1/3) + 2*3**(1/3) * (2*x)**(2/3) - 18*x / (3**(11/6) * (2*x)**(1/3)) - leeDelta, x)
        X = [0, float(mumin[0]), float(mumax[0])]
        Y = [float(mumax[0]), float(mumin[0]), 0]
        return X, Y, mumin, mumax, leeDelta

    resonance = int(resonance[0]) / int(resonance[1])
    X, Y, mumin, mumax, leeDelta = stab(resonance)
    y = np.linspace(0, Y[1], 100)
    x = [float(optimize.brentq(fun, 0, X[2]+.01, args=(y[i], leeDelta))) for i in range(0, len(y))]
    return X, Y, mumin, mumax, x[1:], y[1:]


# # # # # # # # # # # # # # #
# MAKE FIGURE
# # # # # # # # # # # # # # #
def stability_fig_setup(resonance):
    m_lims = {'21': (0, 14.15e-3), '53': (0, 4.15e-3), '32': (0, 4.5e-3), '75': (0, 2.87e-3), '43': (0, 1.42e-3)}
    fig, ax = plt.subplots()
    ax.set_xlim(m_lims[resonance])
    ax.set_ylim(m_lims[resonance])
    ax.set_xlabel('$\mu_1\ [M_1/M_\odot]$')
    ax.set_ylabel('$\mu_2\ [M_2/M_\odot]$')
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 12))
    ax.set_title('{}:{}'.format(*resonance))
    boundary = plot_boundary(resonance, fig, ax)
    return fig, ax, boundary


def plot_boundary(resonance, fig, ax):
    boundary = stability_boundary(resonance)
    limit, = ax.plot(boundary[4], boundary[5], 'b--', label='Analytical limit', linewidth=1)
    ax.plot(boundary[5], boundary[4], 'b--', linewidth=1)
    return boundary


def plot_observed(observed, resonance, fig, ax, boundary, color, label):
    mumin, = boundary[2]
    mumax, = boundary[3]
    for system in observed[resonance]:
        pi_mass = system['pi_mass']
        sig_pi_mass = system['sig_pi_mass']
        po_mass = system['po_mass']
        sig_po_mass = system['sig_po_mass']
        s_mass = system['s_mass']
        sig_s_mass = system['sig_s_mass']
        mu_i, errs_i = mu_error(pi_mass, sig_pi_mass, s_mass, sig_s_mass)
        mu_o, errs_o = mu_error(po_mass, sig_po_mass, s_mass, sig_s_mass)
        ax.errorbar(mu_i, mu_o, xerr=[[errs_i[0]], [errs_i[1]]], yerr=[[errs_o[0]], [errs_o[1]]], ecolor=color, elinewidth=.5, capsize=2, fmt='.', color=color, label=label)
        first_cond = - (mumax - mumin)/mumin*mu_i + mumax
        second_cond = - mumin/(mumax - mumin)*mu_i + mumax*mumin/(mumax - mumin)
        try:
            if mu_o > first_cond and mu_o > second_cond:
                ax.annotate(system['name'], (mu_i, mu_o))
        except:
            pass
