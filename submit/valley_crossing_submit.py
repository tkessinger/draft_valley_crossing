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

s_range = {}
alpha_range = {}
numruns = {}
numsims = {}
sigma_range['161120'] = [0.00001,0.00003,0.0001,0.0003,0.001,0.003,0.01,0.03,0.1] #range of sigma values
rho_range['161120'] = [0] #range of recombination rates
N_range['161120'] = [10000] #range of population size values
delta_range['161120'] = [0.00001,0.00003,0.0001,0.0003,0.001,0.003,0.01,0.03,0.1] #range of deleterious intermediate coefficients
s_range['161120'] = [0.00001,0.00003,0.0001,0.0003,0.001,0.003,0.01,0.03,0.1] #range of selection coefficients for the adaptation
#we will need to think carefully about which ``slices" of parameter space we can explore.
numruns['161120'] = 1
numsamples['161120'] = 10

name = '161120'

sigma_range, rho_range, N_range, delta_range, s_range, numruns, numsamples = sigma_range[name], rho_range[name], N_range[name], delta_range[name], s_range[name], numruns[name], numsamples[name]

if not os.path.exists('output/valley_crossing_sims_'+name):
    os.makedirs('output/valley_crossing_sims_'+name)

shutil.copy('src/simulate_population.py', 'output/valley_crossing_sims_'+name)

for s in s_range:
    for N in N_range:
        for runno in range(numruns):
            call_list = ['sbatch', 'src/valley_crossing.py', '--name', name, '--s', str(s), '--N', str(N), '--numsamples', str(numsamples),
             '--runno', str(runno)]
            call = ' '.join(call_list)
            print call
            #sp.call(call)
            os.system(call)