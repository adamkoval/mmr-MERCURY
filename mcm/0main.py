import os
import argparse
import pandas
import sys
import shutil
from subprocess import Popen

from func import read_envfile

pyenv, bashenv = read_envfile('envfile.txt')

parser = argparse.ArgumentParser()
parser.add_argument('--resonance', '-res',
                    dest='resonance',
                    action='store')
parser.add_argument('--no_runs', '-no',
                    dest='N_runs',
                    action='store')
args = parser.parse_args()

resonance = args.resonance
N_runs = int(args.N_runs)

# Test for nominal resonance directory
if not os.path.exists('../completed'):
    os.mkdir('../completed')

path = '../completed/{}'.format(resonance)
if not os.path.exists(path):
    os.mkdir(path)
    os.mkdir('{}/planets'.format(path))
    os.mkdir('{}/info'.format(path))
    os.mkdir('{}/input'.format(path))
else:
    pass

# Count the number of files which exist already
print('Checking directory structure.')
dirs = ['planets', 'info', 'input']
numbers = {}
for _dir in dirs:
    numbers[_dir] = len(os.listdir('../completed/{}/{}'.format(resonance, _dir)))

if numbers['planets']//2 == numbers['info'] == numbers['input']:
    N_completed = int(numbers['info'])
else:
    print('Please check the numbers of runs in each directory:\n',
          'N_planets: {}'.format(numbers['planets']),
          'N_info: {}'.format(numbers['info']),
          'N_input: {}'.format(numbers['input']))
    sys.exit()

# Main loop
for k in range(N_completed, N_runs+N_completed):
    print('Run = {}'.format(k))

    # Randomizing input
    p_randomize = Popen([pyenv, 'randomize.py', '-res', resonance])
    p_randomize.wait()

    # Cleaning up old files from mercury dir
    p_cleanup = Popen([bashenv, 'cleanup.sh'])
    p_cleanup.wait()

    # Execute and check if Mercury is running
    p_check_mercury = Popen([bashenv, 'check_mercury.sh'])
    p_check_mercury.wait()
    print('Mercury completed, copying files')

    shutil.copyfile('mercury/planet1.aei', '../completed/{}/planets/{}-planet1.aei'.format(resonance, k))
    shutil.copyfile('mercury/planet2.aei', '../completed/{}/planets/{}-planet2.aei'.format(resonance, k))
    shutil.copyfile('mercury/info.out', '../completed/{}/info/{}-info.out'.format(resonance, k))
    shutil.copyfile('mercury/big.in', '../completed/{}/input/{}-big.in'.format(resonance, k))
