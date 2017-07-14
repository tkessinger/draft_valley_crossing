#!/usr/bin/env python
'''
plot_bubbles.py
Author: Taylor Kessinger
Date: June 20, 2017
Description: Glob and plot bubble size distribution.
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

name = '4'

numruns = {}
sigma_range = {}
N_range = {}
delta_range = {}

run_time = {}

if name == '1':

    N_range[name] = [100000]
    sigma_range[name] = [1e-8,1e-6,1e-4,1e-2]
    delta_range[name] = [0.0, 1e-3, 1e-4, 1e-5]
    run_time[name] = 1000000
    numruns[name] = 5
    
    
if name == '2':

    N_range[name] = [1000000]
    sigma_range[name] = [5e-8, 5e-2]
    delta_range[name] = [0.0, 1e-3, 1e-4, 1e-5]
    run_time[name] = 1000000
    numruns[name] = 5

if name == '3':

    N_range[name] = [10000]
    sigma_range[name] = [2e-8, 2e-1]
    delta_range[name] = [0.0, 1e-5, 1e-4, 1e-3, 1e-2]
    run_time[name] = 1000000
    numruns[name] = 5

if name == '4':

    N_range[name] = [10000]
    sigma_range[name] = [2e-8, 2e-1]
    #delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    delta_range[name] = [0.0, 1e-2]
    run_time[name] = 1000000
    numruns[name] = 5

bubble_sizes = {}
bubble_size_sums = {}

for N in N_range[name]:
    mu = 0.1/N
    for sigma in sigma_range[name]:
        for delta in delta_range[name]:
            if not [sigma, delta] in bubble_sizes.keys():
                bubble_sizes[sigma, delta] = []
                bubble_size_sums[sigma, delta] = []

            for runno in range(numruns[name]):
                prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_delta_'+str(delta)+'_sigma_'+str(sigma)
                with open('output/bubble_size_sims_'+name+'/'+prefix+'_weights_'+str(runno)+'.pkl') as weights_file:
                    weights = pickle.load(weights_file)
                    print sigma, delta, runno, len(weights)
                    for wi, weight in enumerate(weights):
                        if weight[-1] != 'bubble closed' and weight[-1] != 'fixed':
                            bubble_sizes[sigma,delta].append(weight)
                            
            for bi, bubble in enumerate(bubble_sizes[sigma,delta]):
                bubble_size_sums[sigma,delta].append(np.sum(bubble))
                
w_theor_bins = np.logspace(1,6,50)
w_2 = w_theor_bins**-2
w_32 = w_theor_bins**-1.5

plt.figure()

logbins = np.logspace(0,6,50)
#for sigma in [sigma_range[name][0],sigma_range[name][-1]]:
#    for delta in [delta_range[name][0],delta_range[name][1]]:



plt.figure()

for sigma in sigma_range[name]:
    #for delta in [delta_range[name][0]]:
    for delta in delta_range[name]:
    
        print sigma, delta, np.average(np.log(bubble_size_sums[sigma,delta])), len(bubble_size_sums[sigma,delta])
        data = np.histogram(bubble_size_sums[sigma,delta],bins=logbins,range=(1,1e6), density=True)
        #data = np.histogram(bubble_size_sums[sigma,delta],bins=20,range=(1,1e6))
        #print data

        plt.plot(data[1][:-1],5*data[0],label='sigma = '+str(sigma) + ', delta = ' + str(delta))
plt.plot(w_theor_bins,w_32,ls='--',c='k')#,label='$w^{-3/2}$')
plt.text(1e4,1e-5,'$w^{-3/2}$')

plt.plot(w_theor_bins,w_2,ls='--',c='k')#,label='$w^{-2}$')
plt.text(1e3,1e-8,'$w^{-2}$')


plt.legend(loc=3)

plt.xlabel('Bubble size $w$')

plt.ylabel('Density $P(w)$')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlim([1e1,1e6])
plt.tight_layout()
plt.show()