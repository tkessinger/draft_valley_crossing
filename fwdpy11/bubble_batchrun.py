#!/usr/bin/env python3
"""
Run bubble_sim.py in batch using SLURM scheduler
"""

import shlex
import subprocess
import itertools as it
import numpy as np

ngens = 50000
burn = 30000
mu = 1e-4
r = 1e-5
reps = 10
seed = 42

pars = ['N', 's']
parvals =  [[1000, 5000, 10000, 50000], list(np.linspace(0.001,0.06,60))]

parsets = list(it.product(*parvals))
npars = len(pars)
nsets = len(parsets)

filename = 'bubble_sim_data.pkl'

for p in range(0,nsets):
    pset = parsets[p]

    filebase = 'bubble_sim'
    simstr = '{}.py '.format(filebase)

    # add parameter values as command line options
    simstr += '--file={} '.format(filename)
    simstr += '--mu={} --r={} --ngens={} --burn={} --reps={} --seed={} '.format(mu,
                                                                                r,
                                                                                ngens,
                                                                                burn,
                                                                                reps,
                                                                                seed)
    simstr += ' '.join(['--{}={}'.format(p[0], p[1]) for p in zip(pars, pset)])

    cmdstr = 'sbatch --job-name=\'bubble_sim\' --output={}.out --wrap=\'{}\''.format(filebase, simstr)

    print(cmdstr)
    cmd = shlex.split(cmdstr)

    proc = subprocess.Popen(cmd)
    proc.wait()
