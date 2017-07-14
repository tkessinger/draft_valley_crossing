#!/usr/bin/env python
'''
plt_asex_crossings.py
Author: Taylor Kessinger
Date: June 16 20, 2017
Description: Glob and plot valley crossing times as a function of input variables.
'''
import sys
sys.path.append('/home/tkessinger/tools/FFPopSim/pkg/python/')
import numpy as np
import FFPopSim as h
import random as rand
import matplotlib.pyplot as plt
from time import time
import argparse
import math
import cPickle as pickle
import glob

import seaborn as sns
sns.set_context('paper', font_scale=1.5)
sns.set_style('white')

name = 'asex_crossings_1'

numruns = {}
sigma_range = {}
N_range = {}
s_range = {}
delta_range = {}

run_time = {}

if name == 'asex_crossings_1':

    N_range[name] = [100000]
    sigma_range[name] = [1e-8,1e-6,1e-4,1e-2]
    delta_range[name] = [0.0, 1e-5, 1e-4,1e-3]
    s_range[name] = [0.01]
    run_time[name] = 1000000
    numruns[name] = 5
    rho = 0.0
    
crossing_times = {}

for N in N_range[name]:
    mu = 0.1/N
    for sigma in sigma_range[name]:
        for s in s_range[name]:
            for delta in delta_range[name]:
                if not [sigma, s, delta] in crossing_times.keys():
                    crossing_times[sigma, s, delta] = []
                
                for runno in range(numruns[name]):
                    prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                    with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                        times = pickle.load(dm_successes_file)
                        for time in times:
                            crossing_times[sigma,s,delta].append(time)
                print sigma, s, delta,len(crossing_times[sigma, s, delta])
                
plt.figure()
for sigma in sigma_range[name]:

    plotted_times = [crossing_times[sigma,s,delta] for delta in delta_range[name]]
    (_, caps, _) = plt.errorbar(delta_range[name],
                                [np.average(x) for x in plotted_times],
                                yerr=[np.std(x) for x in plotted_times], capsize=10,label='$\sigma = $'+str(sigma))

    for cap in caps:
        cap.set_markeredgewidth(1)
plt.legend(loc=2)
plt.xlabel('$\delta$')
plt.ylabel(r"$\tau$")

ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()

plt.show()