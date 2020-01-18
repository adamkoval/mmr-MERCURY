##
#
# Function file for simulation package
#
##
import numpy as np
import os
import sys

def read_envfile(envfile):
    """
    Reads in envfile to assign local environments of
    python and bash.
    In:
        > envfile - (str) path to envfile.txt
    Out:
        > pyenv, bashenv - (str) paths to python and bash
        environments.
    """
    with open(envfile, 'r') as f:
        lines = [line for line in f.readlines() if line[0] is not '#']
        pyenv = lines[0].split()[2]
        bashenv = lines[1].split()[2]
    return pyenv, bashenv


def disp(res_float, fraction, direction=''):
    """
    For the outer planet, sets initial displacement outside of
    or migration distance inside of the resonance under
    consideration.
    In:
        > res_float - (float) fractional value of resonance
        under consideration, i.e., 2/1=2.0, 5/3=1.66667, and
        so on
        > fraction - (float) fraction of distance
        between neighboring resonances by which to
        displace the outer planet
        > direction - (str) 'in' or 'out', depending on
        whether the user sets the initial displacement or
        the final displacement inside of the resonance
    Out:
        > disp - displacement from resonance in T_i/T_o.
    """
    res_floats = [9/7, 4/3, 7/5, 3/2, 5/3, 2/1, 3/1]
    for i in range(1, len(res_floats)-1):
        if res_floats[i] == res_float:
            if direction == 'out':
                disp = fraction * (res_floats[i+1] - res_floats[i])
            elif direction == 'in':
                disp = fraction * (res_floats[i] - res_floats[i-1])
    return disp


def count_completed(res_str):
    """
    Counts number of files present in the 'completed' directory of
    the given resonance.
    In:
        > res_str - (str) specifies the resonance
        under consideration, e.g., '53', '5:3', '5-3'
    Out:
        > N_completed - (int) number of completed runs
        present in the directory.
    """
    print('~~~~~~~~~~~~~~~~~~~~~~~~\n',
          'count_completed():\n',
          'Checking completed directory for the {}:{} resonance.\n'.format(*res_str),
          '~~~~~~~~~~~~~~~~~~~~~~~~\n',)
    dirs = ['planets', 'info', 'input']
    numbers = {}
    for _dir in dirs:
        numbers[_dir] = len(os.listdir('../completed/{}/{}'.format(res_str, _dir)))

    if numbers['planets']//2 == numbers['info'] == numbers['input']:
        N_completed = int(numbers['info'])
        return N_completed
    else:
        print('~~~~~~~~~~~~~~~~~~~~~~~~\n',
              'Please check the numbers of runs in each directory:\n',
              'N_planets: {}\n'.format(numbers['planets']//2),
              'N_info: {}\n'.format(numbers['info']),
              'N_input: {}\n'.format(numbers['input']),
              '~~~~~~~~~~~~~~~~~~~~~~~~\n')
        sys.exit()

def make_directories(res_str):
    def make_if(path):
        if not os.path.exists(path):
            os.mkdir(path)
        else:
            pass

    paths = ['../completed',
             '../completed/{}'.format(res_str),
             '../completed/{}/planets'.format(res_str),
             '../completed/{}/info'.format(res_str),
             '../completed/{}/input'.format(res_str)]

    for path in paths:
        make_if(path)
