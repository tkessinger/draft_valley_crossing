#!/usr/bin/env python
'''
plot_fix_probs.py
Author: Taylor Kessinger
Date: June 20, 2017
Description: Plot fixation probability as a function of s (or delta) and sigma.
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

name = '1'

numruns = {}
sigma_range = {}
N_range = {}
s_range = {}

run_time = {}

if name == '1':

    N_range[name] = [10000]
    sigma_range[name] = [1e-7,1e-6,1e-4,1e-2]
    s_range[name] = [-1e-2,-1e-3,0.0,1e-3,1e-2]
    run_time[name] = 1000000
    numruns[name] = 5
    
if name == '2':

    N_range[name] = [100000]
    sigma_range[name] = [5e-8,5e-2]
    s_range[name] = [-1e-2,-1e-3,0.0,1e-3,1e-2]
    run_time[name] = 1000000
    numruns[name] = 5

if name == '3':

    N_range[name] = [10000]
    sigma_range[name] = [5e-8,5e-2]
    s_range[name] = [-1e-2,-5e-3,-1e-3,0.0,1e-3,5e-3,1e-2]
    run_time[name] = 1000000
    numruns[name] = 10

fix_prob = {}

for N in N_range[name]:
    mu = 0.1/N
    for sigma in sigma_range[name]:
        for s in s_range[name]:
            if not [sigma, s] in fix_prob.keys():
                fix_prob[sigma, s] = []
                failures = 0
                successes = 0
                
            for runno in range(numruns[name]):
                prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_s_'+str(s)+'_sigma_'+str(sigma)
                with open('output/fix_prob_sims_'+name+'/'+prefix+'_weights_'+str(runno)+'.pkl', 'r') as prob_file:
                    probs = pickle.load(prob_file)
                    failures += probs[1]
                    successes += probs[2]
            fix_prob[sigma,s] = 1.0*successes/(successes+failures)
            print sigma, s, successes, failures
                
plt.figure()
for sigma in sigma_range[name]:
    plt.plot(s_range[name], [fix_prob[sigma,s] for s in s_range[name]])
    for s in s_range[name]:
        print sigma, s, fix_prob[sigma,s]
ax=plt.gca()
ax.set_yscale('log')
plt.show()
