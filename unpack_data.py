#
# Utility to unpack all data from old data structure.
# Data structure as:
#   e.g.,   completed/2-1/1/info/1-info.out
#           completed/2-1/1/input/1-big.in
#           completed/2-1/1/planets/1/planet1.aei
#

import os
import argparse
import shutil
import re

parser = argparse.ArgumentParser()
parser.add_argument('--completed_old', '-old',
                    action='store',
                    dest='path_old')
parser.add_argument('--completed_new', '-new',
                    action='store',
                    dest='path_new')
parser.add_argument('--resonance', '-res',
                    action='store',
                    dest='res_str_old')
args = parser.parse_args()

path_old = args.path_old
path_new = args.path_new
resonance =  args.res_str_old
res_str_old = '-'.join((resonance[0], resonance[-1]))
res_str_new = ''.join((resonance[0], resonance[-1]))
three_dirs = ['info', 'input', 'planets']

if not os.path.exists('{}'.format(path_new)):
    os.mkdir(path_new)
if not os.path.exists('{}/{}'.format(path_new, res_str_new)):
    os.mkdir('{}/{}'.format(path_new, res_str_new))
    for _dir in three_dirs:
        os.mkdir('{}/{}/{}'.format(path_new, res_str_new, _dir))

res_dirs = [d for d in os.listdir('{}/{}'.format(path_old, res_str_old)) if re.match('^[0-9]+$', d)]
for res_dir in res_dirs:
    print(res_dir)
    N_old = len(os.listdir('{}/{}/{}/info/'.format(path_old, res_str_old, res_dir)))
    N_new = len(os.listdir('{}/{}/info/'.format(path_new, res_str_new)))
    for n in range(1, N_old+1):
        try:
            shutil.copy('{}/{}/{}/info/{}-info.out'.format(path_old, res_str_old, res_dir, n),
                    '{}/{}/info/{}-info.out'.format(path_new, res_str_new, N_new))
            shutil.copy('{}/{}/{}/input/{}-big.in'.format(path_old, res_str_old, res_dir, n),
                    '{}/{}/input/{}-big.in'.format(path_new, res_str_new, N_new))
            try:
                shutil.copy('{}/{}/{}/planets/{}/planet1.aei'.format(path_old, res_str_old, res_dir, n),
                        '{}/{}/planets/{}-planet1.aei'.format(path_new, res_str_new, N_new))
            except FileNotFoundError:
                open('{}/{}/planets/{}-planet1_empty.aei'.format(path_new, res_str_new, N_new), 'a').close()
            try:
                shutil.copy('{}/{}/{}/planets/{}/planet2.aei'.format(path_old, res_str_old, res_dir, n),
                    '{}/{}/planets/{}-planet2.aei'.format(path_new, res_str_new, N_new))
            except FileNotFoundError:
                open('{}/{}/planets/{}-planet2_empty.aei'.format(path_new, res_str_new, N_new), 'a').close()

            N_new += 1
        except FileNotFoundError:
            pass
