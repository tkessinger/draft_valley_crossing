#!/usr/bin/env python

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-00:30:00

# Set name of job shown in squeue
#SBATCH --job-name src/branching_cluster.py

# Request CPU resources
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
'''
valley_crossing_submit.py
Author: Taylor Kessinger
Date: November 19, 2016
Description: Submit script for valley crossing simulations.
'''

import os
import shutil
import sys
import numpy as np

name = sys.argv[1]

numruns = {}
sigma_range = {}
N_range = {}
s_range = {}

run_time = {}

if name == '1':

    N_range[name] = [10000]
    sigma_range[name] = [1e-7,1e-6,1e-4,1e-2]
    s_range[name] = [-1e-2,-1e-3,0,1e-3,1e-2]
    run_time[name] = 1000000
    numruns[name] = 5

if name == '2':

    N_range[name] = [100000]
    sigma_range[name] = [5e-8,5e-2]
    s_range[name] = [-1e-2,-1e-3,0,1e-3,1e-2]
    run_time[name] = 1000000
    numruns[name] = 5

if name == '3':

    N_range[name] = [10000]
    sigma_range[name] = [5e-8,5e-2]
    s_range[name] = [-1e-2,-5e-3,-1e-3,0.0,1e-3,5e-3,1e-2]
    run_time[name] = 1000000
    numruns[name] = 10


if not os.path.exists('output/fix_prob_sims_'+name):
    os.makedirs('output/fix_prob_sims_'+name)

shutil.copy('src/fixation_probability.py', 'output/fix_prob_sims_'+name)

for N in N_range[name]:
    for sigma in sigma_range[name]:
        for s in s_range[name]:
            for runno in range(numruns[name]):
                call_list = ['sbatch', 'src/fixation_probability.py', '--name', name, '--N', str(N), '--mu', str(0.1/N),
                             '--sigma', str(sigma), '--s', str(s), '--runtime', str(run_time[name]), '--runno', str(runno)]
                
                call = ' '.join(call_list)
                print call
                #sp.call(call)
                os.system(call)
