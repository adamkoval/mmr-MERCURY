import os
import argparse
import pandas
import sys
import shutil
from subprocess import Popen

import func as fn

pyenv, bashenv = fn.read_envfile('envfile.txt')

parser = argparse.ArgumentParser()
parser.add_argument('--resonance', '-res',
                    dest='resonance',
                    action='store')
parser.add_argument('--no_runs', '-no',
                    dest='N_runs',
                    action='store')
parser.add_argument('--process_no', '-pno',
                    dest='process_no',
                    action='store')
args = parser.parse_args()

res_str = args.resonance
N_runs = int(args.N_runs)
pno = str(args.process_no)

# Test for nominal resonance directory
fn.make_directories(res_str)

# Count the number of files which exist already
N_completed = fn.count_completed(res_str)

# Main loop
k = 0
while k < N_runs:
#for k in range(N_completed, N_runs+N_completed):
    print('Run = {}'.format(k))

    # Randomizing input
    p_randomize = Popen([pyenv, 'randomize.py', '-res', res_str, '-pno', pno])
    p_randomize.wait()

    # Cleaning up old files from mercury dir
    p_cleanup = Popen([bashenv, 'cleanup.sh', pno])
    p_cleanup.wait()

    # Execute and check if Mercury is running
    p_check_mercury = Popen([bashenv, 'check_mercury.sh', pno])
    p_check_mercury.wait()
    N_completed = fn.count_completed(res_str)
    print('N_completed = {}'.format(N_completed))
    print('Mercury completed, copying files')

    shutil.copyfile('mercury_{}/planet1.aei'.format(pno), '../completed/{}/planets/{}-planet1.aei'.format(res_str, N_completed))
    shutil.copyfile('mercury_{}/planet2.aei'.format(pno), '../completed/{}/planets/{}-planet2.aei'.format(res_str, N_completed))
    shutil.copyfile('mercury_{}/info.out'.format(pno), '../completed/{}/info/{}-info.out'.format(res_str, N_completed))
    shutil.copyfile('mercury_{}/big.in'.format(pno), '../completed/{}/input/{}-big.in'.format(res_str, N_completed))
    k += 1
