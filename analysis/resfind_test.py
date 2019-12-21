import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import resfind_funcs as resfind
import general_funcs as gen

"""
Key:
    [1]  pl_hostname - host star
    [2]  pl_letter - planet letter
    [3]  pl_pnum - number of planets
    [4]  pl_orbper - orbital period [days]
    [5]  pl_orbpererr1 - max error [days]
    [6]  pl_orbpererr2 - min error [days]
    [7]  pl_orbperlim - "Orbital period Limit flag."
    [8]  pl_bmassj - mass (in Mass or M*sin(i))
    [9]  pl_bmassjerr1 - max error
    [10]  pl_bmassjerr2 - min error
    [11] pl_bmassjlim - "Planet Mass or M*sin(i) Limit flag."
    [12] pl_bmassprov - mass provenance
    [13] st_mass - stellar mass [Msol]
    [14] st_masserr1 - max error
    [15] st_masserr2 - min error
"""

planets = pd.read_csv('/home/akoval/uni/MMR_PH50007/old/summer/codee/data_planets/planets.csv', skiprows=23, header=0, delim_whitespace=False)
tol = 15

#pmasses, smasses, names, mtypes, merrors = resfind.find_resonances(planets, tol, output=False)

res_floats = [2/1, 5/3, 3/2, 7/5, 4/3]
res_strs = ['21', '53', '32', '75', '43']

"""
def find_resonances(planets, tol, output):

"""
g = 0
systems = {}
planets_headers = [key for key in planets]
for i in range(1, len(planets)):
    if planets['pl_hostname'][i] == planets['pl_hostname'][g]:
        continue
    else:
        sys_name = planets['pl_hostname'][g]
        systems[sys_name] = planets[g:i]
        g = i

#for arr in ['pmasses', 'smasses', 'names', 'mtypes', 'merrors']:
#    globals()[arr] = {'21':[], '53':[], '32':[], '75':[], '43':[]}
#sys_data = {'names': names, 'pmasses': pmasses, 'smasses': smasses,
#            'mtypes': mtypes, 'merrors': merrors}

sys_data = {'21': [], '53': [], '32': [], '75': [], '43': []}

sys_names = [key for key in systems]

for i in range(len(systems)):
    curr_sys = systems[sys_names[i]]
    curr_sys_info = {'name': [], 'pimass': [], 'pimasserrs': [],
                     'pomass': [], 'pomasserrs': [], 'smass': [],
                     'smasserrs': [], 'mtype': []}
    for j in range(len(curr_sys)-1):
            for k in range(j+1, len(curr_sys)):
                massj = np.asarray(curr_sys)[j]
                massk = np.asarray(curr_sys)[k]
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

                    """
                    def check_tol(tol, pmasses, smasses, names,
                                  mtypes, merrors, systems, j, k, i):
                    """
                    p_o = np.asarray(curr_sys['pl_orbper'])[outer]
                    p_o_maxerr = np.asarray(curr_sys['pl_orbpererr1'])[outer]
                    p_o_minerr = np.asarray(curr_sys['pl_orbpererr2'])[outer]
                    p_i = np.asarray(curr_sys['pl_orbper'])[inner]
                    p_i_maxerr = np.asarray(curr_sys['pl_orbpererr1'])[inner]
                    p_i_minerr = np.asarray(curr_sys['pl_orbpererr2'])[inner]

                    errs = []
                    for resonance in res_floats:
                        """
                        def tolerance(resonance, p_o, min_err_o, max_err_o,
                                    p_i, min_err_i, max_err_i):
                        """
                        max_period = [p_o+p_o_maxerr, p_i+p_i_maxerr]
                        min_period = [p_o+p_o_minerr, p_i+p_i_minerr]
                        max_ratio = max_period[0] / min_period[1]
                        min_ratio = min_period[0] / max_period[1]
                        max_per_diff = max_ratio / resonance
                        min_per_diff = min_ratio / resonance
                        upp_per_err = abs(1 - max_per_diff)
                        low_per_err = abs(1 - min_per_diff)

                        if upp_per_err > low_per_err:
                            errs.append(upp_per_err * 100)
                        else:
                            errs.append(low_per_err * 100)
                        """
                        end tolerance()
                        """
                    if min(errs) < tol:
                        #print((max_ratio+min_ratio)/2)
                        for l in range(len(res_floats)):
                            if min(errs) == errs[l]:
                                """
                                def append_list(sys_data[i], systems, iteration, lesser, larger):
                                """
                                curr_sys_info = {'name': [], 'pimass': [], 'pimasserrs': [],
                                                     'pomass': [], 'pomasserrs': [], 'smass': [],
                                                     'smasserrs': [], 'mtype': []}

                                curr_sys_info['name'] = '{} {}, {}'.format(
                                        np.asarray(curr_sys['pl_hostname'])[inner],
                                        np.asarray(curr_sys['pl_letter'])[inner],
                                        np.asarray(curr_sys['pl_letter'])[outer])
                                curr_sys_info['pimass'] = np.asarray(curr_sys['pl_bmassj'])[inner]
                                curr_sys_info['pimasserrs'] = (
                                        np.asarray(curr_sys['pl_bmassjerr2'])[inner],
                                        np.asarray(curr_sys['pl_bmassjerr1'])[inner])
                                curr_sys_info['pomass'] = np.asarray(curr_sys['pl_bmassj'])[outer]
                                curr_sys_info['pomasserrs'] = (
                                        np.asarray(curr_sys['pl_bmassjerr2'])[outer],
                                        np.asarray(curr_sys['pl_bmassjerr1'])[outer])
                                curr_sys_info['smass'] = np.asarray(curr_sys['st_mass'])[inner]
                                curr_sys_info['smasserrs'] = (
                                        np.asarray(curr_sys['st_masserr2'])[inner],
                                        np.asarray(curr_sys['st_masserr1'])[inner])
                                curr_sys_info['mtype'] = np.asarray(curr_sys['pl_bmassprov'])[inner]
                                curr_sys_info['piper'] = np.asarray(curr_sys['pl_orbper'])[inner]
                                curr_sys_info['pipererrs'] = (
                                        np.asarray(curr_sys['pl_orbpererr2'])[inner],
                                        np.asarray(curr_sys['pl_orbpererr1'])[inner])
                                curr_sys_info['poper'] = np.asarray(curr_sys['pl_orbper'])[outer]
                                curr_sys_info['popererrs'] = (
                                        np.asarray(curr_sys['pl_orbpererr2'])[outer],
                                        np.asarray(curr_sys['pl_orbpererr1'])[outer])

                                sys_data[res_strs[l]].append(curr_sys_info)
                                """
                                end append_list()
                                """
                            else:
                                pass
                    """
                    end check_tol()
                    """

# Need to test for resonances.
plotting_data = {}
for resonance in res_strs:
    for system in sys_data[resonance]:
        p_ratio = system['poper']/system['piper']
        p_ratio_min = (system['poper'] + system['popererrs'][0]) / (system['piper'] + system['pipererrs'][1])
        p_ratio_max = (system['poper'] + system['popererrs'][1]) / (system['piper'] + system['pipererrs'][0])
        plotting_data[system['name']] = (p_ratio, p_ratio_min, p_ratio_max)

y = [i for i in range(len(plotting_data))]
y_labels = [key for key in plotting_data]

fig, ax = plt.subplots()
for i, resonance in enumerate(res_floats):
    ax.axvline(resonance)
    ax.text(resonance, 1, res_strs[i])
for i, system in enumerate(plotting_data):
    x, xmin, xmax = plotting_data[system]
    #plt.errorbar(x, y[i], xerr=([xmin, xmax]))
    plt.errorbar(x, y[i], marker='*')
plt.show()

for resonance in res_strs:
    fig, ax = plt.subplots()
    for system in sys_data[resonance]:
        ypos = system['pomass']
        xpos = system['pimass']
        ax.plot(xpos, ypos, '*')
    ax.set_title(resonance)
    plt.show()
