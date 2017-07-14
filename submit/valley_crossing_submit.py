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

s_range = {}
numruns = {}
numsims = {}
sigma_range = {}
N_range = {}
rho_range = {}
delta_range = {}
mu_range = {}

numsamples = {}

sigma_range['161120'] = [0.00001] #range of sigma values
rho_range['161120'] = [0] #range of recombination rates
N_range['161120'] = [1000] #range of population size values
delta_range['161120'] = [0.00001] #range of deleterious intermediate coefficients
s_range['161120'] = [0.00001] #range of selection coefficients for the adaptation
#we will need to think carefully about which ``slices" of parameter space we can explore.
numruns['161120'] = 1
numsamples['161120'] = 1

sigma_range['weissman_2009_1'] = [1e-10] #range of sigma values
rho_range['weissman_2009_1'] = [0] #range of recombination rates
N_range['weissman_2009_1'] = [10**2,2*10**2, 3*10**2,5*10**2,10**3,2*10**3, 3*10**3, 5*10**3, 10**4, 2*10**4, 3*10**4, 5*10**4, 10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
delta_range['weissman_2009_1'] = [2*10**-4] #range of deleterious intermediate coefficients
s_range['weissman_2009_1'] = [0.1] #range of selection coefficients for the adaptation
mu_range['weissman_2009_1'] = [10**-5]
numruns['weissman_2009_1'] = 10
numsamples['weissman_2009_1'] = 10

sigma_range['weissman_2009_2'] = [1e-10] #range of sigma values
rho_range['weissman_2009_2'] = [0] #range of recombination rates
N_range['weissman_2009_2'] = [10**2,2*10**2, 3*10**2,5*10**2,10**3,3*10**3,5*10**3,10**4,3*10**4,5*10**4,10**5, 3*10**5,5*10**5, 10**6, 3*10**6, 5*10**6,10**7]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
delta_range['weissman_2009_2'] = [7*10**-3] #range of deleterious intermediate coefficients
s_range['weissman_2009_2'] = [0.1] #range of selection coefficients for the adaptation
mu_range['weissman_2009_2'] = [10**-5]
numruns['weissman_2009_2'] = 10
numsamples['weissman_2009_2'] = 10

sigma_range['weissman_2009_3'] = [1e-10] #range of sigma values
rho_range['weissman_2009_3'] = [0] #range of recombination rates
N_range['weissman_2009_3'] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
delta_range['weissman_2009_3'] = [-0.01,-0.005,-.001,0.0,.001,.005,.01] #range of deleterious intermediate coefficients
s_range['weissman_2009_3'] = [0.1] #range of selection coefficients for the adaptation
mu_range['weissman_2009_3'] = [10**-6]#[4*10**-5]
numruns['weissman_2009_3'] = 10
numsamples['weissman_2009_3'] = 10

sigma_range['weissman_2009_4'] = [1e-10] #range of sigma values
rho_range['weissman_2009_4'] = [0] #range of recombination rates
N_range['weissman_2009_4'] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
delta_range['weissman_2009_4'] = [3*10**-3] #range of deleterious intermediate coefficients
s_range['weissman_2009_4'] = [0.1] #range of selection coefficients for the adaptation
mu_range['weissman_2009_4'] = [10**-5, 3*10**-5, 10**-4, 3*10**-4, 10**-3, 3*10**-3, 10**-2]
numruns['weissman_2009_4'] = 10
numsamples['weissman_2009_4'] = 10

sigma_range['weissman_2010_1'] = [1e-10] #range of sigma values
rho_range['weissman_2010_1'] = 2.0*np.array([10**-6,10**-5,3*10**-5,10**-4,3*10**-4,10**-3,3*10**-3,10**-2,3*10**-2,10**-1,3*10**-1]) #range of recombination rates
N_range['weissman_2010_1'] = [10**5] #range of population size values
delta_range['weissman_2010_1'] = [10**-5] #range of deleterious intermediate coefficients
s_range['weissman_2010_1'] = [0.05] #range of selection coefficients for the adaptation
mu_range['weissman_2010_1'] = [5*10**-7]
numruns['weissman_2010_1'] = 5
numsamples['weissman_2010_1'] = 10

sigma_range['weissman_2010_2'] = [1e-10] #range of sigma values
rho_range['weissman_2010_2'] = 2.0*np.array([10**-6,10**-5,3*10**-5,10**-4,3*10**-4,10**-3,3*10**-3,10**-2,3*10**-2,10**-1,3*10**-1]) #range of recombination rates
N_range['weissman_2010_2'] = [10**8] #range of population size values
delta_range['weissman_2010_2'] = [5*10**-5] #range of deleterious intermediate coefficients
s_range['weissman_2010_2'] = [0.05] #range of selection coefficients for the adaptation
mu_range['weissman_2010_2'] = [2*10**-9]
numruns['weissman_2010_2'] = 5
numsamples['weissman_2010_2'] = 10

sigma_range['weissman_2010_3'] = [1e-10] #range of sigma values
rho_range['weissman_2010_3'] = 2.0*np.array([10**-6,10**-5,3*10**-5,10**-4,3*10**-4,10**-3,3*10**-3,10**-2,3*10**-2,10**-1,3*10**-1]) #range of recombination rates
N_range['weissman_2010_3'] = [10**6] #range of population size values
delta_range['weissman_2010_3'] = [5*10**-5] #range of deleterious intermediate coefficients
s_range['weissman_2010_3'] = [0.05] #range of selection coefficients for the adaptation
mu_range['weissman_2010_3'] = [2*10**-9]
numruns['weissman_2010_3'] = 5
numsamples['weissman_2010_3'] = 10

sigma_range['weissman_2010_4'] = [1e-10] #range of sigma values
rho_range['weissman_2010_4'] = 2.0*np.array([10**-6,10**-5,3*10**-5,10**-4,3*10**-4,10**-3,3*10**-3,10**-2,3*10**-2,10**-1,3*10**-1]) #range of recombination rates
N_range['weissman_2010_4'] = [10**6] #range of population size values
delta_range['weissman_2010_4'] = [5*10**-5] #range of deleterious intermediate coefficients
s_range['weissman_2010_4'] = [0.05] #range of selection coefficients for the adaptation
mu_range['weissman_2010_4'] = [2*10**-7]
numruns['weissman_2010_4'] = 5
numsamples['weissman_2010_4'] = 10


sigma_range, rho_range, N_range, delta_range, s_range, mu_range, numruns, numsamples = sigma_range[name], rho_range[name], N_range[name], delta_range[name], s_range[name], mu_range[name], numruns[name], numsamples[name]

if not os.path.exists('output/valley_crossing_sims_'+name):
    os.makedirs('output/valley_crossing_sims_'+name)

shutil.copy('src/simulate_population.py', 'output/valley_crossing_sims_'+name)

for s in s_range:
    for N in N_range:
        for sigma in sigma_range:
            for rho in rho_range:
                for mu in mu_range:
                    for delta in delta_range:
                        for runno in range(numruns):
                            call_list = ['sbatch', 'src/simulate_population.py', '--name', name, '--s', str(s), '--N', str(N), '--sigma', str(sigma), '--rho', str(rho),
                                         '--mu', str(mu), '--delta', str(delta), '--numsamples', str(numsamples),
                                         '--runno', str(runno)]
                            call = ' '.join(call_list)
                            print call
                            #sp.call(call)
                            os.system(call)