#!/usr/bin/env python
'''
plt_asex_crossings.py
Author: Taylor Kessinger
Date: June 16, 2017
Description: Glob and plot valley crossing times as a function of input variables.
'''
import sys
sys.path.append('/home/tkessinger/tools/FFPopSim/pkg/python/')
sys.path.append('/home/tkessinger/Documents/draft_valley_crossing/FFPopSim/pkg/python/')
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
sns.set_style('ticks')

#recently used plots, possibly paper worthy:
#in all cases N*mu ~= 0.3
#var_sigma_6: variable sigma, N = 10**7
#var_sigma_7: variable sigma, N = 10**6
#var_sigma_8: variable sigma, N = 10**5


name = 'var_sigma_sequential_fix_N_4_delta_-5'

numruns = {}
sigma_range = {}
N_range = {}
s_range = {}
delta_range = {}
mu_range = {}
rho_range = {}

run_time = {}

if name == 'bands':
    graphtype = 'bands'
    fixed_s = False

if name == 'asex_crossings_1':

    N_range[name] = [100000]
    sigma_range[name] = [1e-8,1e-6,1e-4,1e-2]
    delta_range[name] = [0.0, 1e-5, 1e-4,1e-3]
    s_range[name] = [0.01]
    run_time[name] = 1000000
    numruns[name] = 5
    rho = 0.0
    graphtype = 'sigma_delta'
    
if name == 'asex_crossings_2':

    N_range[name] = [100000]
    sigma_range[name] = [1e-8,1e-6,1e-4,1e-2]
    delta_range[name] = [0.0, 1e-3]
    s_range[name] = [0.01]
    run_time[name] = 1000000
    numruns[name] = 5
    rho_range[name] = [0]
    fixed_s = False
    graphtype = 'sigma_delta'

if name == 'asex_crossings_2_fixed_s':

    N_range[name] = [100000]
    sigma_range[name] = [5e-8,5e-6,5e-4,5e-2]
    delta_range[name] = [0.0, 1e-3]
    s_range[name] = [0.01]
    run_time[name] = 1000000
    numruns[name] = 5
    rho_range[name] = [0]
    fixed_s = True
    graphtype = 'sigma_delta'

if name == 'asex_weissman_2009_1':
    N_range[name] = [10**2,3*10**2,10**3,3*10**3, 10**4, 3*10**4, 10**5, 3*10**5, 10**6]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [1e-10] #range of sigma values
    delta_range[name] = [2*10**-4] #range of deleterious intermediate coefficients
    s_range[name] = [0.1] #range of selection coefficients for the adaptation
    mu_range[name] = [10**-5]
    numruns[name] = 5
    run_time[name] = 1000000
    rho_range[name] = [0]
    rho = 0.0
    fixed_s = False
    graphtype = 'N'
    
if name == 'asex_weissman_2009_2':
    N_range[name] = [10**2,3*10**2,10**3,3*10**3, 10**4, 3*10**4, 10**5, 3*10**5, 10**6]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [1e-10] #range of sigma values
    delta_range[name] = [7*10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [0.1] #range of selection coefficients for the adaptation
    mu_range[name] = [10**-5]
    numruns[name] = 5
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'N'

if name == 'asex_weissman_2009_3':

    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [1e-10] #range of sigma values
    delta_range[name] = [-0.01,-0.005,-.002,-.001,0.0,.001,.002,.005,.01] #range of deleterious intermediate coefficients #range of deleterious intermediate coefficients
    s_range[name] = [0.1] #range of selection coefficients for the adaptation
    mu_range[name] = [10**-6]
    numruns[name] = 5
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'delta'

if name == 'sex_weissman_2010_1':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [1e-10] #range of sigma values
    delta_range[name] = [10**-5] #range of deleterious intermediate coefficients
    s_range[name] = [0.05] #range of selection coefficients for the adaptation
    mu_range[name] = [5*10**-7]
    numruns[name] = 5
    rho_range[name] = 2*np.array([10**-6, 3*10**-6, 10**-5, 3*10**-5, 10**-4, 3*10**-4, 10**-3, 3*10**-3, 10**-2, 3*10**-2, 10**-1, 3*10**-1])
    run_time[name] = 1000000
    fixed_s = False
    graphtype = 'rho'


if name == 'sex_weissman_2010_2':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [1e-10] #range of sigma values
    delta_range[name] = [10**-5] #range of deleterious intermediate coefficients
    s_range[name] = [0.05] #range of selection coefficients for the adaptation
    mu_range[name] = [5*10**-7]
    numruns[name] = 5
    rho_range[name] = 2*np.array([10**-6, 3*10**-6, 10**-5, 3*10**-5, 10**-4, 3*10**-4, 10**-3, 3*10**-3, 10**-2, 3*10**-2, 10**-1, 3*10**-1])
    run_time[name] = 1000000
    fixed_s = False
    graphtyoe = 'rho'

if name == 'var_sigma_1':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    delta_range[name] = [7*10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [0.1] #range of selection coefficients for the adaptation
    mu_range[name] = [10**-5]
    numruns[name] = 5
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'
    
if name == 'var_sigma_2':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])

    delta_range[name] = [7*10**-2] #range of deleterious intermediate coefficients
    s_range[name] = [0.1] #range of selection coefficients for the adaptation
    mu_range[name] = [10**-5]
    numruns[name] = 5
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

# if name == 'var_sigma_3':
#     N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
#     sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
#     delta_range[name] = [7*10**-2] #range of deleterious intermediate coefficients
#     s_range[name] = [0.1] #range of selection coefficients for the adaptation
#     mu_range[name] = [10**-6]
#     numruns[name] = 5
#     run_time[name] = 1000000
#     rho = 0.0
#     fixed_s = False
    
if name == 'var_sigma_3':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    delta_range[name] = [7*10**-2] #range of deleterious intermediate coefficients
    s_range[name] = [0.1] #range of selection coefficients for the adaptation
    mu_range[name] = [3*10**-6]
    numruns[name] = 5
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

if name == 'var_sigma_4':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    delta_range[name] = [7*10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [0.01] #range of selection coefficients for the adaptation
    mu_range[name] = [3*10**-6]
    numruns[name] = 5
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'
    
if name == 'var_sigma_5':
    N_range[name] = [5*10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    delta_range[name] = [7*10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [0.01] #range of selection coefficients for the adaptation
    mu_range[name] = [3*10**-6]
    numruns[name] = 5
    run_time[name] = 10000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'
    
if name == 'sigma_delta_1':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [5e-6, 5e-2] #range of sigma values
    delta_range[name] = [10**-5, 10**-4,10**-3,10**-2] #range of deleterious intermediate coefficients
    s_range[name] = [10**-4, 10**-3,10**-2,10**-1] #range of selection coefficients for the adaptation
    mu_range[name] = [3*10**-6]
    numruns[name] = 5
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 's_delta'
    
if name == 'sigma_delta_2':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [5e-6, 5e-2] #range of sigma values
    delta_range[name] = [10**-5, 10**-4,10**-3,10**-2] #range of deleterious intermediate coefficients
    #s_range[name] = [10**-4, 10**-3,10**-2,10**-1] #range of selection coefficients for the adaptation
    s_range[name] = [10**-1] #range of selection coefficients for the adaptation
    mu_range[name] = [3*10**-6]
    numruns[name] = 5
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 's_delta'

    
if name == 'var_sigma_6':
    N_range[name] = [10**7]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [10**-2] #range of selection coefficients for the adaptation
    mu_range[name] = [3*10**-8]
    numruns[name] = 5
    run_time[name] = 10000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

if name == 'var_sigma_7':
    N_range[name] = [10**6]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [10**-2] #range of selection coefficients for the adaptation
    mu_range[name] = [3*10**-7]
    numruns[name] = 5
    run_time[name] = 10000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

    
if name == 'var_sigma_8':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    #sigma_range[name] = [10**x for x in range(-7,-1)] + [3*10**x for x in range(-7,-1)] #range of sigma values
    sigma_range[name] = [10**x for x in range(-5,-1)] + [3*10**x for x in range(-5,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [10**-2] #range of selection coefficients for the adaptation
    mu_range[name] = [3*10**-6]
    numruns[name] = 5
    run_time[name] = 10000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

if name == 'var_sigma_sequential_fix_N_4_delta_-5':
    N_range[name] = [1*10**4]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    sigma_range[name] = np.array(sigma_range[name])
    delta_range[name] = [10**-5] #range of deleterious intermediate coefficients
    s_range[name] = [10**-2] #range of selection coefficients for the adaptation
    mu_range[name] = [1*10**-5]
    numruns[name] = 3
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

if name == 'var_sigma_sequential_fix_N_3_delta_-5':
    N_range[name] = [1*10**3]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    sigma_range[name] = np.array(sigma_range[name])
    delta_range[name] = [10**-5] #range of deleterious intermediate coefficients
    s_range[name] = [10**-2] #range of selection coefficients for the adaptation
    mu_range[name] = [1*10**-5]
    numruns[name] = 3
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'
    
if name == 'var_sigma_deterministic_tunnel':
    N_range[name] = [1*10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    sigma_range[name] = np.array(sigma_range[name])
    delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [10**-2] #range of selection coefficients for the adaptation
    mu_range[name] = [1*10**-4]
    numruns[name] = 6
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

if name == 'var_sigma_deterministic_fix':
    N_range[name] = [1*10**6]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
    sigma_range[name] = np.sort(sigma_range[name])
    sigma_range[name] = np.array(sigma_range[name])
    delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [10**-2] #range of selection coefficients for the adaptation
    mu_range[name] = [1*10**-3]
    numruns[name] = 6
    run_time[name] = 1000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

if name == 'var_sigma_tunneling':
    N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)]    
    sigma_range[name] = np.sort(sigma_range[name])
    delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
    s_range[name] = [10**-2] #range of selection coefficients for the adaptation
    mu_range[name] = [1*10**-6]
    numruns[name] = 3
    run_time[name] = 10000000
    rho = 0.0
    fixed_s = False
    graphtype = 'sigma'

crossing_times = {}
if fixed_s:
    inferred_sigmas = {}

inferred_Uvals = {}

def plt_sigma_delta():
    for N in N_range[name]:
        mu = 0.1/N
        for sigma in sigma_range[name]:
            for s in s_range[name]:
                for delta in delta_range[name]:
    
                    if not [sigma, s, delta] in crossing_times.keys():
                        crossing_times[sigma, s, delta] = []
                    
                    for runno in range(numruns[name]):
                        if fixed_s:
                            prefix = 'N_'+str(N)+'_smean_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                        else:
                            prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                        
                        with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                            times = pickle.load(dm_successes_file)
                            for time in times:
                                crossing_times[sigma,s,delta].append(time)
                        if fixed_s:
                            if not sigma in inferred_sigmas.keys():
                                inferred_sigmas[sigma] = []
                            with open('output/valley_crossing_sims_'+name+'/'+prefix+'_sigma_hist_'+str(runno)+'.pkl', 'r') as sigma_hist_file:
                                svals = pickle.load(sigma_hist_file)
                                for sval in svals:
                                    inferred_sigmas[sigma].append(sval)
                    print sigma, s, delta,len(crossing_times[sigma, s, delta])
    if fixed_s:
        print sigma_range[name], [np.average(np.array(inferred_sigmas[sigma])[:,0]) for sigma in sigma_range[name]]
    plt.figure()
    for sigma in sigma_range[name]:
      
        plotted_times = [crossing_times[sigma,s,delta] for delta in delta_range[name]]
        if fixed_s:
            sigma = np.average(np.array(inferred_sigmas[sigma])[:,0])
        (_, caps, _) = plt.errorbar(np.array(delta_range[name])+1e-6,
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
    
def plt_N_dependence():
    for N in N_range[name]:
        for mu in mu_range[name]:
            for sigma in sigma_range[name]:
                for s in s_range[name]:
                    for delta in delta_range[name]:
                        params = (N,mu,sigma,s,delta)
        
                        if not params in crossing_times.keys():
                            crossing_times[params] = []
                        
                        for runno in range(numruns[name]):
                            if fixed_s:
                                prefix = 'N_'+str(N)+'_smean_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            else:
                                prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            
                            with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                times = pickle.load(dm_successes_file)
                                for time in times:
                                    crossing_times[params].append(time)
                            if fixed_s:
                                if not sigma in inferred_sigmas.keys():
                                    inferred_sigmas[sigma] = []
                                with open('output/valley_crossing_sims_'+name+'/'+prefix+'_sigma_hist_'+str(runno)+'.pkl', 'r') as sigma_hist_file:
                                    svals = pickle.load(sigma_hist_file)
                                    for sval in svals:
                                        inferred_sigmas[sigma].append(sval)
                        print N, sigma, s, delta,len(crossing_times[params])
        if fixed_s:
            print sigma_range[name], [np.average(np.array(inferred_sigmas[sigma])[:,0]) for sigma in sigma_range[name]]
    plt.figure()
    for sigma in sigma_range[name]:
      
        plotted_times = [crossing_times[N,mu,sigma,s,delta] for N in N_range[name]]
        if fixed_s:
            sigma = np.average(np.array(inferred_sigmas[sigma])[:,0])
        (_, caps, _) = plt.errorbar(np.array(N_range[name]),
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

def plt_delta_dependence():
    for N in N_range[name]:
        for mu in mu_range[name]:
            for sigma in sigma_range[name]:
                for s in s_range[name]:
                    for delta in delta_range[name]:
                        params = (N,mu,sigma,s,delta)
        
                        if not params in crossing_times.keys():
                            crossing_times[params] = []
                        
                        for runno in range(numruns[name]):
                            if fixed_s:
                                prefix = 'N_'+str(N)+'_smean_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            else:
                                prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            
                            with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                times = pickle.load(dm_successes_file)
                                for time in times:
                                    crossing_times[params].append(time)
                            if fixed_s:
                                if not sigma in inferred_sigmas.keys():
                                    inferred_sigmas[sigma] = []
                                with open('output/valley_crossing_sims_'+name+'/'+prefix+'_sigma_hist_'+str(runno)+'.pkl', 'r') as sigma_hist_file:
                                    svals = pickle.load(sigma_hist_file)
                                    for sval in svals:
                                        inferred_sigmas[sigma].append(sval)
                        print N, sigma, s, delta,len(crossing_times[params])
        if fixed_s:
            print sigma_range[name], [np.average(np.array(inferred_sigmas[sigma])[:,0]) for sigma in sigma_range[name]]
    plt.figure()
    for sigma in sigma_range[name]:
      
        plotted_times = [crossing_times[N,mu,sigma,s,delta] for delta in delta_range[name]]
        if fixed_s:
            sigma = np.average(np.array(inferred_sigmas[sigma])[:,0])
        (_, caps, _) = plt.errorbar(np.array(delta_range[name]),
                                    [np.average(x) for x in plotted_times],
                                    yerr=[np.std(x) for x in plotted_times], capsize=10,label='$\sigma = $'+str(sigma))
      
        for cap in caps:
            cap.set_markeredgewidth(1)
    plt.legend(loc=2)
    plt.xlabel('$\delta$')
    plt.ylabel(r"$\tau$")
      
    ax=plt.gca()
    #ax.set_xscale('log')
    ax.set_yscale('log')
    plt.tight_layout()
      
    plt.show()
    
def plt_rho_dependence():
    for N in N_range[name]:
        for mu in mu_range[name]:
            for sigma in sigma_range[name]:
                for s in s_range[name]:
                    for delta in delta_range[name]:
                        for rho in rho_range[name]:
                            params = (N,mu,sigma,s,delta, rho)
            
                            if not params in crossing_times.keys():
                                crossing_times[params] = []
                            
                            for runno in range(numruns[name]):
                                if fixed_s:
                                    prefix = 'N_'+str(N)+'_smean_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                                else:
                                    prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                                
                                with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                    times = pickle.load(dm_successes_file)
                                    for time in times:
                                        crossing_times[params].append(time)
                                if fixed_s:
                                    if not sigma in inferred_sigmas.keys():
                                        inferred_sigmas[sigma] = []
                                    with open('output/valley_crossing_sims_'+name+'/'+prefix+'_sigma_hist_'+str(runno)+'.pkl', 'r') as sigma_hist_file:
                                        svals = pickle.load(sigma_hist_file)
                                        for sval in svals:
                                            inferred_sigmas[sigma].append(sval)
                            print N, sigma, s, delta,len(crossing_times[params])
        if fixed_s:
            print sigma_range[name], [np.average(np.array(inferred_sigmas[sigma])[:,0]) for sigma in sigma_range[name]]
    plt.figure()
    for sigma in sigma_range[name]:
      
        plotted_times = [crossing_times[N,mu,sigma,s,delta,rho] for rho in rho_range[name]]
        if fixed_s:
            sigma = np.average(np.array(inferred_sigmas[sigma])[:,0])
        (_, caps, _) = plt.errorbar(np.array(rho_range[name]),
                                    2*np.array([np.average(x) for x in plotted_times]),
                                    yerr=[np.std(x) for x in plotted_times], capsize=10,label='$\sigma = $'+str(sigma))
      
        for cap in caps:
            cap.set_markeredgewidth(1)
    plt.legend(loc=2)
    plt.xlabel(r'$\rho$')
    plt.ylabel(r"$\tau$")
      
    ax=plt.gca()
    ax.set_xscale('log')
    #ax.set_yscale('log')
    plt.tight_layout()
      
    plt.show()

def plt_sigma_dependence():
    for N in N_range[name]:
        for mu in mu_range[name]:
            for sigma in sigma_range[name]:
                for s in s_range[name]:
                    for delta in delta_range[name]:
                        params = (N,mu,sigma,s,delta)
        
                        if not params in crossing_times.keys():
                            crossing_times[params] = []
                        
                        for runno in range(numruns[name]):
                            if fixed_s:
                                prefix = 'N_'+str(N)+'_smean_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            else:
                                prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            
                            with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                times = pickle.load(dm_successes_file)
                                for time in times:
                                    crossing_times[params].append(time)
                            if fixed_s:
                                if not sigma in inferred_sigmas.keys():
                                    inferred_sigmas[sigma] = []
                                with open('output/valley_crossing_sims_'+name+'/'+prefix+'_sigma_hist_'+str(runno)+'.pkl', 'r') as sigma_hist_file:
                                    svals = pickle.load(sigma_hist_file)
                                    for sval in svals:
                                        inferred_sigmas[sigma].append(sval)
                        print N, sigma, s, delta,len(crossing_times[params])
        if fixed_s:
            print sigma_range[name], [np.average(np.array(inferred_sigmas[sigma])[:,0]) for sigma in sigma_range[name]]
    plt.figure()
      
    plotted_times = [crossing_times[N,mu,sigma,s,delta] for sigma in sigma_range[name]]
    if fixed_s:
        sigma = np.average(np.array(inferred_sigmas[sigma])[:,0])
    (_, caps, _) = plt.errorbar(np.array(sigma_range[name]),
                                [np.average(x) for x in plotted_times],
                                yerr=[np.std(x) for x in plotted_times], capsize=10)
      
    for cap in caps:
        cap.set_markeredgewidth(1)
    plt.legend(loc=2)
    plt.xlabel('$\sigma$')
    plt.ylabel(r"$\tau$")
      
    ax=plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.tight_layout()
      
    plt.show()
    
    plt.figure()
    plotted_times = [np.log10(crossing_times[N,mu,sigma,s,delta]) for sigma in sigma_range[name]]
    #plt.boxplot(plotted_times, positions=np.log10(np.array(sigma_range[name])), manage_xticks=False,widths=0.10,meanline=True)
    plt.boxplot(plotted_times, positions=np.log10(np.array(sigma_range[name])), manage_xticks=False,widths=0.08,meanline=True)
    
    plt.legend(loc=2)
    plt.xlabel(r'$\log_{10} (\sigma)$')
    plt.ylabel(r"$\log_{10} (\tau)$")
    #plt.xticks([-5,-4,-3,-2,-1])
      
    ax=plt.gca()
    plt.tight_layout()
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.set_xticks(ax.get_xticks()[::4]+ax.get_xticks()[-1])
    #plt.tight_layout()
    plt.show()


def plt_s_delta_dependence():
    for N in N_range[name]:
        for mu in mu_range[name]:
            for sigma in sigma_range[name]:
                for s in s_range[name]:
                    for delta in delta_range[name]:
                        params = (N,mu,sigma,s,delta)
        
                        if not params in crossing_times.keys():
                            crossing_times[params] = []
                        
                        for runno in range(numruns[name]):
                            if fixed_s:
                                prefix = 'N_'+str(N)+'_smean_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            else:
                                prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            
                            with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                times = pickle.load(dm_successes_file)
                                for time in times:
                                    crossing_times[params].append(time)
                            if fixed_s:
                                if not sigma in inferred_sigmas.keys():
                                    inferred_sigmas[sigma] = []
                                with open('output/valley_crossing_sims_'+name+'/'+prefix+'_sigma_hist_'+str(runno)+'.pkl', 'r') as sigma_hist_file:
                                    svals = pickle.load(sigma_hist_file)
                                    for sval in svals:
                                        inferred_sigmas[sigma].append(sval)
                            if not params in inferred_Uvals.keys():
                                inferred_Uvals[params] = []
                            tmp_Uvals = []
                            with open('output/valley_crossing_sims_'+name+'/'+prefix+'_effective_mut_rate_'+str(runno)+'.pkl', 'r') as Uval_hist_file:
                                Uvals = pickle.load(Uval_hist_file)
                                for Uval in Uvals:
                                    tmp_Uvals.append(Uval)
                        inferred_Uvals[params] = 1.0*np.sum(tmp_Uvals)/N/200
                        print N, sigma, s, delta,len(crossing_times[params])
    for sigma in sigma_range[name]:
        for s in [s_range[name][0],s_range[name][-1]]:
            plotted_times = [crossing_times[N,mu,sigma,s,delta] for delta in delta_range[name]]
            (_, caps, _) = plt.errorbar(np.array(delta_range[name]),
                            [np.average(x) for x in plotted_times],
                            yerr=[np.std(x) for x in plotted_times], capsize=10, label=(sigma,s))
  
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
    
    plt.figure()
    colors = [['red','pink'], ['blue', 'cyan']]
    ci = 0
    for sigma in sigma_range[name]:
        plotted_times = [np.log10(crossing_times[N,mu,sigma,s,delta]) for delta in delta_range[name]]
        #plt.boxplot(plotted_times, positions=np.log10(np.array(sigma_range[name])), manage_xticks=False,widths=0.10,meanline=True)
        bp = plt.boxplot(plotted_times, positions=np.log10(np.array(delta_range[name])), manage_xticks=False,widths=0.15,meanline=True)#,patch_artist=True)#, label=r'$\sigma = '+str(sigma))
        
        for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bp[element], color=colors[ci][0])
    
        for patch in bp['boxes']:
            patch.set(color = colors[ci][1]) 
        ci += 1
    
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='red', lw=1),
                Line2D([0], [0], color='blue', lw=1)]

    plt.legend(custom_lines, [r'$\sigma = 5 \times 10^{-6}$',r'$\sigma = 5 \times 10^{-2}$'],loc=4)
    
    #plt.legend(loc=2)
    plt.xlabel(r'$\log_{10} (\delta)$')
    plt.ylabel(r"$\log_{10} (\tau)$")
    #plt.xticks([-5,-4,-3,-2,-1])
      
    ax=plt.gca()
    plt.tight_layout()
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.set_xticks(ax.get_xticks()[::4]+ax.get_xticks()[-1])
    #plt.tight_layout()
    plt.show()

def plt_bands_indv(line_label,name):
    if name == 'bands_1_1':
        N_range[name] = [10**4/2]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = 2.0*np.array(sigma_range[name])
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-5]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False
    
    if name == 'bands_1_2':
        N_range[name] = [10**4]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = np.array(sigma_range[name])
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-5]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False
    
    if name == 'bands_1_3':
        N_range[name] = [2*10**4]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = np.array(sigma_range[name])/2.0
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-5]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False
    
    if name == 'bands_2_1':
        N_range[name] = [10**5/2]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = 2.0*np.array(sigma_range[name])
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-6]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False
    
    if name == 'bands_2_2':
        N_range[name] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = np.array(sigma_range[name])
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-6]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False
    
    if name == 'bands_2_3':
        N_range[name] = [2*10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = np.array(sigma_range[name])/2.0
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-6]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False
    
    if name == 'bands_3_1':
        N_range[name] = [10**6/2]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = 2.0*np.array(sigma_range[name])
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-6]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False
    
    if name == 'bands_3_2':
        N_range[name] = [10**6]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = np.array(sigma_range[name])
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-7]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False
    
    if name == 'bands_3_3':
        N_range[name] = [2*10**6]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
        sigma_range[name] = [10**x for x in range(-5,-1)] + [2*10**x for x in range(-5,-1)] + [8*10**x for x in range(-5,-1)] + [5*10**x for x in range(-5,-1)] #range of sigma values
        sigma_range[name] = np.sort(sigma_range[name])
        sigma_range[name] = np.array(sigma_range[name])/2.0
        delta_range[name] = [10**-3] #range of deleterious intermediate coefficients
        s_range[name] = [10**-2] #range of selection coefficients for the adaptation
        mu_range[name] = [3*10**-7]
        numruns[name] = 3
        run_time[name] = 1000000
        rho = 0.0
        fixed_s = False

    for N in N_range[name]:
        for mu in mu_range[name]:
            for sigma in sigma_range[name]:
                for s in s_range[name]:
                    for delta in delta_range[name]:
                        params = (N,mu,sigma,s,delta)
        
                        if not params in crossing_times.keys():
                            crossing_times[params] = []
                        
                        for runno in range(numruns[name]):
                            if fixed_s:
                                prefix = 'N_'+str(N)+'_smean_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            else:
                                prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            
                            with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                times = pickle.load(dm_successes_file)
                                for time in times:
                                    crossing_times[params].append(time)
                        print name, N, sigma, s, delta,len(crossing_times[params])
      
    plotted_times = [crossing_times[N,mu,sigma,s,delta] for sigma in sigma_range[name]]
    [x.sort() for x in plotted_times]
    A = [len(x) for x in plotted_times]
    print A
    asymmetric_error = [[-np.average(x) + x[np.int(len(x)/4)] for x in plotted_times],[-np.average(x) + x[np.int(len(x)/4)] for x in plotted_times]]
    (_, caps, _) = plt.errorbar(N*np.array(sigma_range[name]),
                                [np.average(x) for x in plotted_times],
                                #yerr=[np.std(x) for x in plotted_times], capsize=10, label=line_label)
                                yerr = asymmetric_error,capsize=10,label=line_label)
    
    plt.tight_layout()
    plt.show()
    for cap in caps:
        cap.set_markeredgewidth(0.5)    

def plt_bands():
    plt.figure()
    name = 'bands_1_1'
    plt_bands_indv('low N',name)
    name = 'bands_1_2'
    plt_bands_indv('medium N',name)
    name = 'bands_1_3'
    plt_bands_indv('high N',name)
    plt.legend(loc=2)
    plt.xlabel('$N\sigma$')
    plt.ylabel(r"$\tau$")
    plt.title('$N = 10^4$')
      
    ax=plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.tight_layout()

    plt.figure()
    name = 'bands_2_1'
    plt_bands_indv('low N',name)
    name = 'bands_2_2'
    plt_bands_indv('medium N',name)
    name = 'bands_2_3'
    plt_bands_indv('high N',name)
    plt.legend(loc=2)
    plt.xlabel('N$\sigma$')
    plt.ylabel(r"$\tau$")
    plt.title('$N = 10^5$')

    
    ax=plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.tight_layout()
    
    plt.figure()
    name = 'bands_3_1'
    plt_bands_indv('low N',name)
    name = 'bands_3_2'
    plt_bands_indv('medium N',name)
    name = 'bands_3_3'
    plt_bands_indv('high N',name)
    plt.legend(loc=2)
    plt.xlabel('N$\sigma$')
    plt.ylabel(r"$\tau$")
    plt.title('$N = 10^6$')

    
    ax=plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.tight_layout()
    
    plt.show()

#       
#     ax=plt.gca()
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     plt.tight_layout()
#       
#     plt.show()       


if graphtype == 'sigma':
    plt_sigma_dependence()
elif graphtype == 's_delta':
    plt_s_delta_dependence()
elif graphtype == 'rho':
    plt_rho_dependence()
elif graphtype == 'sigma_delta':
    plt_sigma_delta()
elif graphtype == 'N':
    plt_N_dependence()
elif graphtype == 'delta':
    plt_delta_dependence()
elif graphtype == 'bands':
    plt_bands()