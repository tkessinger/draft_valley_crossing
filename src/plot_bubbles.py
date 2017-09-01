#!/usr/bin/env python
'''
plot_bubbles.py
Author: Taylor Kessinger
Date: June 20, 2017
Description: Glob and plot bubble size distribution.
'''
import sys
sys.path.append('/home/tkessinger/tools/FFPopSim/pkg/python/')
sys.path.append('/home/tkessinger/.local/lib/python3.5/site-packages')
import numpy as np
#import FFPopSim as h
import random as rand
import matplotlib.pyplot as plt
from time import time
import argparse
import math
import cPickle as pickle
import glob
import powerlaw as pl

import seaborn as sns
sns.set_context('paper', font_scale=1.5)
sns.set_style('white')

name = '10'

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
    
if name == '5':

    N_range[name] = [10000]
    sigma_range[name] = [2e-8, 2e-1]
    #delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    delta_range[name] = [0.0, 1e-2]
    run_time[name] = 50000
    numruns[name] = 10

if name == '6':

    s_mean_range = {}
    N_range[name] = [10000]
    sigma_range[name] = [1e-6, 1e-2]
    #delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    delta_range[name] = [0.0, 1e-2]
    run_time[name] = 50000
    numruns[name] = 10
    
if name == '7':

    s_mean_range = {}
    N_range[name] = [50000]
    sigma_range[name] = [1e-6, 1e-2]
    #delta_range[name] = [0, 1e-3, 1e-4, 1e-5]
    delta_range[name] = [0.0, 1e-2]
    run_time[name] = 200000
    numruns[name] = 10
    
if name == '8':

    N_range[name] = [100000]
    sigma_range[name] = [2e-8, 2e-1]
    delta_range[name] = [0.0, 1e-2]
    run_time[name] = 1000000
    numruns[name] = 10

if name == '9':

    N_range[name] = [100000]
    sigma_range[name] = [2e-8, 2e-1]
    delta_range[name] = [0.0, 1e-2]
    run_time[name] = 1000000
    numruns[name] = 10
    
if name == '10':

    N_range[name] = [100000]
    sigma_range[name] = [2e-8, 2e-1]
    delta_range[name] = [0.0, 1e-2]
    run_time[name] = 1000000
    numruns[name] = 10

bubble_sizes = {}
bubble_size_sums = {}
bubble_durations = {}

for N in N_range[name]:
    mu = 0.1/N
    for sigma in sigma_range[name]:
        for delta in delta_range[name]:
            if not (sigma, delta) in bubble_sizes.keys():
                bubble_sizes[sigma, delta] = []
                bubble_size_sums[sigma, delta] = []
                bubble_durations[sigma,delta] = []

            for runno in range(numruns[name]):
                prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_delta_'+str(delta)+'_sigma_'+str(sigma)
                with open('output/bubble_size_sims_'+name+'/'+prefix+'_weights_'+str(runno)+'.pkl', 'rb') as weights_file:
                    #print(weights_file)
                    weights = pickle.load(weights_file)
                    new_weights = []
                    [new_weights.extend(x) for x in weights]
#    weights = [y for y in x for x in weights]
                    weights = new_weights
                    print(sigma, delta, runno, len(weights))
                    for wi, weight in enumerate(weights):
                        if weight[-1] != 'bubble closed' and weight[-1] != 'fixed':
                            bubble_sizes[sigma,delta].append(weight)
                            bubble_durations[sigma,delta].append(len(weight))
                            
            for bi, bubble in enumerate(bubble_sizes[sigma,delta]):
                bubble_size_sums[sigma,delta].append(np.sum(bubble))
            print(sigma, delta, np.average(bubble_size_sums[sigma,delta]),np.median(bubble_size_sums[sigma,delta]), np.average(bubble_durations[sigma,delta]))
                
w_theor_bins = np.logspace(2,7,50)
w_2 = w_theor_bins**-2
w_32 = w_theor_bins**-1.5
w_1 = w_theor_bins**-1
w_12 = w_theor_bins**-0.5


logbins = np.logspace(2,7,50)
#for sigma in [sigma_range[name][0],sigma_range[name][-1]]:
#    for delta in [delta_range[name][0],delta_range[name][1]]:



plt.figure()

for sigma in sigma_range[name]:
    #for delta in [delta_range[name][0]]:
    for delta in delta_range[name]:
    
        print(sigma, delta, np.average(np.log(bubble_size_sums[sigma,delta])), len(bubble_size_sums[sigma,delta]))
        data = np.histogram(bubble_size_sums[sigma,delta],bins=logbins,range=(10,1e6))
        #data = np.histogram(bubble_size_sums[sigma,delta])#,range=(1,1e5))
        hist_dat = 1.0*data[0]
        #hist_dat*=w_32[0]/hist_dat[0]
        hist_dat /= (data[1][1:] - data[1][:-1])
        plt.plot(data[1][:-1],hist_dat,label='$\sigma$ = '+str(sigma) + ', $\delta$ = '+str(delta))

        #data = np.histogram(bubble_size_sums[sigma,delta],bins=20,range=(1,1e6))
        #print data

        #plt.plot(data[1][:-1],5*data[0],label='sigma = '+str(sigma) + ', delta = ' + str(delta))

plt.legend(loc=3)

plt.xlabel('Bubble size $w$')

plt.ylabel('Density $P(w)$')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlim([1e2,1e7])
plt.tight_layout()
plt.show()



plt.figure()

for sigma in sigma_range[name]:
    #for delta in [delta_range[name][0]]:
    for delta in delta_range[name]:
    
        print(sigma, delta, np.average(np.log(bubble_size_sums[sigma,delta])), len(bubble_size_sums[sigma,delta]))
        data = np.histogram(bubble_size_sums[sigma,delta],bins=logbins)#,range=(10,1e6))
        #data = np.histogram(bubble_size_sums[sigma,delta])#,range=(1,1e5))
        hist_dat = 1.0*data[0]
        #hist_dat /= (data[1][1:] - data[1][:-1])
        cumhist = np.cumsum(hist_dat)[-1] - np.cumsum(hist_dat)
        #cumhist /= (data[1][1:] - data[1][:-1])
        cumhist*=w_1[0]/cumhist[0]
        
        plt.plot(data[1][:-1],cumhist,label='$s$ = '+str(sigma) + ', $\delta$ = '+str(delta))

        #data = np.histogram(bubble_size_sums[sigma,delta],bins=20,range=(1,1e6))
        #print data

        #plt.plot(data[1][:-1],5*data[0],label='sigma = '+str(sigma) + ', delta = ' + str(delta))
plt.plot(w_theor_bins,w_1,ls='--',c='k')#,label='$w^{-3/2}$')
#plt.text(1e4,1e-5,'$w^{-3/2}$')
plt.plot(w_theor_bins,w_12,ls='--',c='k')#,label='$w^{-2}$')


plt.legend(loc=3)

plt.xlabel('Bubble size $w$')

plt.ylabel('cdf of $w$')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlim([1e2,1e7])
plt.tight_layout()
plt.show()

# 
# plt.figure()
# for delta in delta_range[name]:
#     data_1 = np.histogram(bubble_size_sums[2e-8,delta],bins=logbins)
#     data_2 = np.histogram(bubble_size_sums[2e-1,delta],bins=logbins)
# 
#     hist_dat = 1.0*data_1[0]
#     cumhist_1 = np.cumsum(hist_dat)[-1] - np.cumsum(hist_dat)
# 
#     hist_dat = 1.0*data_2[0]
#     cumhist_2 = np.cumsum(hist_dat)[-1] - np.cumsum(hist_dat)
#     
#     plt.plot(data[1][:-1],data_1[0] - data_2[0], label=delta)
# plt.show()
# ax=plt.gca()
# ax.set_xscale('log')
# plt.legend(loc=3)
# 
# 
# plt.figure()
# for delta in delta_range[name]:
#     data_1 = np.histogram(bubble_size_sums[2e-8,delta],bins=logbins)
#     data_2 = np.histogram(bubble_size_sums[2e-1,delta],bins=logbins)
# 
#     hist_dat = 1.0*data_1[0]
#     cumhist_1 = np.cumsum(hist_dat)[-1] - np.cumsum(hist_dat)
# 
#     hist_dat = 1.0*data_2[0]
#     cumhist_2 = np.cumsum(hist_dat)[-1] - np.cumsum(hist_dat)
#     
#     plt.plot(data[1][:-1],cumhist_1-cumhist_2, label=delta)
# plt.show()
# ax=plt.gca()
# ax.set_xscale('log')
# plt.legend(loc=3)


def plot_weight_powerlaw(weights, fit = None, loc = 3, labelstr = ''):
    if fit == None:
        fit = pl.Fit(weights, discrete=True)

    fit.plot_pdf(color = 'b', label = 'empirical '+ labelstr)
    fit.power_law.plot_pdf(color = 'b', linestyle = '--',
                               label=r'power law fit: $\alpha$ = ' + str(np.round(fit.alpha, 3)))
    plt.legend(loc=3)

    return fit

plt.figure()
for sigma in sigma_range[name]:
    #for delta in [delta_range[name][0]]:
    for delta in delta_range[name]:
    
        print(sigma, delta, np.average(np.log(bubble_size_sums[sigma,delta])), len(bubble_size_sums[sigma,delta]))
        data = bubble_size_sums[sigma,delta]
        #data = np.histogram(bubble_size_sums[sigma,delta])#,range=(1,1e5))
        plot_weight_powerlaw(data, labelstr = str((sigma, delta)))
    
    
plt.show()



'''
plt.figure()

for sigma in sigma_range[name]:
    #for delta in [delta_range[name][0]]:
    for delta in delta_range[name]:
    
        data = np.histogram(bubble_durations[sigma,delta],bins=logbins)#,range=(10,1e6))
        #data = np.histogram(bubble_size_sums[sigma,delta])#,range=(1,1e5))
        hist_dat = 1.0*data[0]
        #hist_dat /= (data[1][1:] - data[1][:-1])
        cumhist = np.cumsum(hist_dat)[-1] - np.cumsum(hist_dat)
        #cumhist /= (data[1][1:] - data[1][:-1])
        cumhist*=w_1[0]/cumhist[0]
        
        plt.plot(data[1][:-1],cumhist,label='$\sigma$ = '+str(sigma) + ', $\delta$ = '+str(delta))

        #data = np.histogram(bubble_size_sums[sigma,delta],bins=20,range=(1,1e6))
        #print data

        #plt.plot(data[1][:-1],5*data[0],label='sigma = '+str(sigma) + ', delta = ' + str(delta))
plt.plot(w_theor_bins,w_1,ls='--',c='k')#,label='$w^{-3/2}$')
#plt.text(1e4,1e-5,'$w^{-3/2}$')

plt.plot(w_theor_bins,w_12,ls='--',c='k')#,label='$w^{-2}$')
#plt.text(1e3,1e-8,'$w^{-2}$')
plt.plot(w_theor_bins,w_32,ls='--',c='k')#,label='$w^{-3/2}$')



plt.legend(loc=3)

plt.xlabel('Bubble size $w$')

plt.ylabel('Density $P(w)$')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlim([1e1,1e6])
plt.tight_layout()
plt.show()
'''