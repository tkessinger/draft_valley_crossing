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
bubble_sizes_submit.py
Author: Taylor Kessinger
Date: June 15, 2017
Description: Submit script for bubble size simulations.
'''

import os
import shutil
import sys
import numpy as np

name = sys.argv[1]

numruns = {}
sigma_range = {}
N_range = {}
delta_range = {}

run_time = {}

if name == '1':

    N_range[name] = [100000]
    sigma_range[name] = [1e-8,1e-6,1e-4,1e-2]
    delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    run_time[name] = 1000000
    numruns[name] = 5

if name == '2':

    N_range[name] = [1000000]
    sigma_range[name] = [5e-8, 5e-2]
    delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    run_time[name] = 1000000
    numruns[name] = 5
    
if name == '3':

    N_range[name] = [10000]
    sigma_range[name] = [2e-8, 2e-1]
    #delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    delta_range[name] = [1e-2]
    run_time[name] = 1000000
    numruns[name] = 5
    
if name == '4':

    N_range[name] = [10000]
    sigma_range[name] = [2e-8, 2e-1]
    #delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    delta_range[name] = [0.0, 1e-2]
    run_time[name] = 1000000
    numruns[name] = 5


if not os.path.exists('output/bubble_size_sims_'+name):
    os.makedirs('output/bubble_size_sims_'+name)

shutil.copy('src/bubble_sizes.py', 'output/bubble_size_sims_'+name)

for N in N_range[name]:
    for sigma in sigma_range[name]:
        for delta in delta_range[name]:
            for runno in range(numruns[name]):
                call_list = ['sbatch', 'src/bubble_sizes.py', '--name', name, '--N', str(N), '--mu', str(0.1/N),
                             '--sigma', str(sigma), '--delta', str(delta), '--runtime', str(run_time[name]), '--runno', str(runno)]
                
                call = ' '.join(call_list)
                print call
                #sp.call(call)
                os.system(call)
