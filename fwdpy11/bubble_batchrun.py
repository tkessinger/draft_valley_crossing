#!/usr/bin/env python3
"""
Run bubble_sim.py in batch using SLURM scheduler
"""

import shlex
import subprocess
import itertools as it
import numpy as np

ngens = 10000
theta = 0.4
r = 1e-6
reps = 10
seed = 42

pars = ['N', 's']
parvals = [[1000, 5000, 10000, 50000], list(np.linspace(0.002, 0.1, 50))]

parsets = list(it.product(*parvals))
npars = len(pars)
nsets = len(parsets)

outputbase = 'bubble_sim_data_18-03-2018'

for p in range(0,nsets):
    pset = dict(zip(pars, parsets[p]))
    
    simstr = 'bubble_sim.py '
    # add parameter values as command line options
    simstr += '--file={} '.format(outputbase + '.pkl')
    simstr += (
        '--N={} --s={} --mu={} --r={} --ngens={} --burn={} --reps={} --seed={} '
        .format(pset['N'],
                pset['s'],
                theta/(4*pset['N']),
                r,
                ngens,
                5*pset['N'],
                reps,
                seed))

    cmdstr = (
        'sbatch --job-name=\'bubble_sim\' --output={}_{}.out --wrap=\'{}\''
        .format(outputbase,
                str(p+1).zfill(len(str(nsets))),
                simstr))

    print(cmdstr)
    # cmd = shlex.split(cmdstr)

    # proc = subprocess.Popen(cmd)
    # proc.wait()
