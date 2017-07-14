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
asex_crossings_submit.py
Author: Taylor Kessinger
Date: June 10, 2017
Description: Submit valley crossing simulations for fixed sigma.
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
delta_range = {}

run_time = {}

if name == 'asex_crossings_1':

    N_range[name] = [100000]
    sigma_range[name] = [1e-8,1e-6,1e-4,1e-2]
    delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    s_range[name] = [0.01]
    run_time[name] = 1000000
    numruns[name] = 5

if not os.path.exists('output/valley_crossing_sims_'+name):
    os.makedirs('output/valley_crossing_sims_'+name)

shutil.copy('src/simulate_population.py', 'output/valley_crossing_sims_'+name)

for N in N_range[name]:
    for sigma in sigma_range[name]:
        for s in s_range[name]:
            for delta in delta_range[name]:
                for runno in range(numruns[name]):
                    call_list = ['sbatch', 'src/simulate_population.py', '--name', name, '--N', str(N), '--mu', str(0.1/N),
                                 '--sigma', str(sigma), '--delta', str(delta), '--s', str(s), '--runtime', str(run_time[name]), '--runno', str(runno)]
                    
                    call = ' '.join(call_list)
                    print call
                    #sp.call(call)
                    os.system(call)
