##
#
# Function file for simulation package
#
##
import numpy as np

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
