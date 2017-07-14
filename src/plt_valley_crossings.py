#!/usr/bin/env python
'''
simulate_population.py
Author: Taylor Kessinger
Date: February 20, 2017
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

def weissman_asexual_times(N, s, delta, mu):
    N, s, delta, mu = np.array(N), np.array(s), np.array(delta), np.array(mu)
    times = []
    for N in N_range:
        for s in s_range:
            for delta in delta_range:
                for mu in mu_range:
                    mu = 2.0*mu
                    #p_1 = (-delta + np.sqrt(delta**2 + 4*mu*s))/2
                    #times.append(1.0/N/mu/p_1 + 1.0/N/mu/s + 0.577/s)
                    if N > 2*s/(np.pi*mu**2) and N > delta/mu**2: #deterministic fixation
                        times.append(np.log((s+delta)/N*mu**2))
                        print N, s, delta, mu
                        print "deterministic"
                    elif N > 1/mu and N > 2*delta**2/(np.pi*mu**2*s): #neutral semi-deterministic tunneling
                        times.append(np.sqrt(np.pi/(2*N*mu**2*s))+1.0/s*(np.log(1-s)-np.log(s)))
                        print N, s, delta, mu
                        print "neutral semi-deterministic"
                    elif N > 1.0/np.sqrt(mu*s) and delta < 2*np.sqrt(mu*s): #neutral stochastic tunneling
                        #times.append(1.0/(N*mu*np.sqrt(mu*s)))
                        p1 = (-delta+np.sqrt(delta**2+4*mu*s))/2
                        times.append(1.0/(N*mu*p1))
                        print N, s, delta, mu
                        print "neutral tunneling"
                    elif N > 1.0/delta*np.log(delta**2/(mu*s)): #deleterious tunneling
                        #times.append(delta/(N*mu**2*s))
                        p1 = (-delta+np.sqrt(delta**2+4*mu*s))/2
                        times.append(1.0/(N*mu*p1))
                        print N, s, delta, mu
                        print "deleterious tunneling"
                    elif N > 1.0/delta: #deleterious sequential fixation
                        rho = (np.exp(delta)-1)/(np.exp(N*delta)-1)
                        times.append(1.0/(N*mu*rho))
                        print N, s, delta, mu
                        print "deleterious sequential fixation"
                    else: #neutral sequential fixation, lower left quadrant
                        times.append(1.0/mu)
                        print N, s, delta, mu
                        print "neutral sequential fixation"
                        
                        
                        #times.append(delta/(N*mu**2*s))
    return times

#for name in ['weissman_2009_1', 'weissman_2009_2', 'weissman_2009_3', 'weissman_2009_4']:
for name in ['weissman_2010_1', 'weissman_2010_4']:#, 'weissman_2009_3', 'weissman_2009_4']:
#for name in ['weissman_2009_4']:

    print name

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
    N_range['weissman_2009_2'] = [10**3,3*10**3,10**4,3*10**4,10**5, 3*10**5, 10**6, 3*10**6, 10**7]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    delta_range['weissman_2009_2'] = [7*10**-3] #range of deleterious intermediate coefficients
    s_range['weissman_2009_2'] = [0.1] #range of selection coefficients for the adaptation
    #mu_range['weissman_2009_2'] = [10**-6]
    mu_range['weissman_2009_2'] = [10**-5]
    numruns['weissman_2009_2'] = 10
    numsamples['weissman_2009_2'] = 10
    
    
    sigma_range['weissman_2009_3'] = [1e-10] #range of sigma values
    rho_range['weissman_2009_3'] = [0] #range of recombination rates
    N_range['weissman_2009_3'] = [10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
    delta_range['weissman_2009_3'] = [-0.01,-0.005,-.001,.001,.005,.01] #range of deleterious intermediate coefficients
    s_range['weissman_2009_3'] = [0.1] #range of selection coefficients for the adaptation
    mu_range['weissman_2009_3'] = [1e-6]#[4*10**-5]
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
    numruns['weissman_2010_1'] = 10
    numsamples['weissman_2010_1'] = 10
    
    sigma_range['weissman_2010_2'] = [1e-10] #range of sigma values
    rho_range['weissman_2010_2'] = 2.0*np.array([10**-6,10**-5,3*10**-5,10**-4,3*10**-4,10**-3,3*10**-3,10**-2,3*10**-2,10**-1,3*10**-1]) #range of recombination rates
    N_range['weissman_2010_2'] = [10**8] #range of population size values
    delta_range['weissman_2010_2'] = [5*10**-5] #range of deleterious intermediate coefficients
    s_range['weissman_2010_2'] = [0.05] #range of selection coefficients for the adaptation
    mu_range['weissman_2010_2'] = [2*10**-9]
    numruns['weissman_2010_2'] = 10
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

# filemask = 'output/valley_crossing_sims_'+name+'/*'+'.pkl'
# 
# files = glob.glob(filemask)
# 
# for fname in files:


    valley_crossing_times = {}
    
    s_range, N_range, sigma_range, rho_range, mu_range, delta_range, numruns = s_range[name], N_range[name], sigma_range[name], rho_range[name], mu_range[name], delta_range[name], numruns[name]
    
    
    for s in s_range:
        for N in N_range:
            for sigma in sigma_range:
                for rho in rho_range:
                    for mu in mu_range:
                        for delta in delta_range:
                            if not [N,sigma,mu,s,delta,rho] in valley_crossing_times.keys():
                                valley_crossing_times[N,sigma,mu,s,delta,rho] = []
                            prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)

#                            prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_s_'+str(s)+'_delta_'+str(delta)
                            for runno in range(numruns):
                                with open('output/valley_crossing_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                    crossing_times = pickle.load(dm_successes_file)
                                    for vt, vtime in enumerate(crossing_times):
                                        if vt == 0:
                                            valley_crossing_times[N,sigma,mu,s,delta,rho].append(crossing_times[0])
                                        else:
                                            valley_crossing_times[N,sigma,mu,s,delta,rho].append(crossing_times[vt])
    
    if name == 'weissman_2009_1':
        plt.figure()
        plotted_times = [valley_crossing_times[N,sigma,mu,s,delta,rho] for N in N_range]
        (_, caps, _) = plt.errorbar(N_range, [np.average(x) for x in plotted_times], yerr=[np.std(x) for x in plotted_times], capsize=10)

        for cap in caps:
                cap.set_markeredgewidth(1)
        theor = weissman_asexual_times(np.array(N_range), s, delta, mu)
        plt.plot(N_range, theor, ls = '--')
        plt.xlabel('$N$')
        plt.ylabel(r'$\tau$')
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()
        plt.show()
        
    if name == 'weissman_2009_2':
        plt.figure()
        plotted_times = [valley_crossing_times[N,sigma,mu,s,delta,rho] for N in N_range]
        (_, caps, _) = plt.errorbar(N_range, [np.average(x) for x in plotted_times], yerr=[np.std(x) for x in plotted_times], capsize=10)

        for cap in caps:
                cap.set_markeredgewidth(1)
        theor = weissman_asexual_times(np.array(N_range), s, delta, mu)
        plt.plot(N_range, theor, ls = '--')        
        plt.xlabel('$N$')
        plt.ylabel(r'$\tau$')
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()

        plt.show()
        
    if name == 'weissman_2009_3':
        plt.figure()
        plotted_times = [valley_crossing_times[N,sigma,mu,s,delta,rho] for delta in delta_range]
        (_, caps, _) = plt.errorbar(delta_range, [np.average(x) for x in plotted_times], yerr=[np.std(x) for x in plotted_times], capsize=10)

        for cap in caps:
                cap.set_markeredgewidth(1)
        print [np.average(x) for x in plotted_times]
        theor = weissman_asexual_times(N, s, np.array(delta_range), mu)
        plt.plot(delta_range, theor, ls = '--')
        plt.xlabel('$\delta$')
        plt.ylabel(r'$\tau$')
        ax=plt.gca()
        
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        plt.tight_layout()

        plt.show()
        
    if name == 'weissman_2009_4':
        plt.figure()
        plotted_times = [valley_crossing_times[N,sigma,mu,s,delta,rho] for mu in mu_range]
        (_, caps, _) = plt.errorbar(mu_range, [np.average(x) for x in plotted_times], yerr=[np.std(x) for x in plotted_times], capsize=10)

        for cap in caps:
                cap.set_markeredgewidth(1)
        theor = weissman_asexual_times(N, s, delta, mu_range)
        plt.plot(mu_range, theor, ls = '--')
        plt.xlabel('$\mu$')
        plt.ylabel(r"$\tau$")
        
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()

        plt.show()

    if name == 'weissman_2010_1':
        plt.figure()
        plotted_times = [valley_crossing_times[N,sigma,mu,s,delta,rho] for rho in rho_range]
        plt.errorbar(rho_range, [np.average(x) for x in plotted_times], yerr=[np.std(x) for x in plotted_times])#, capsize=10)

        #theor = weissman_asexual_times(N, s, delta, mu_range)
        #plt.plot(rho_range, theor, ls = '--')
        #plt.xlabel('$\rho$')
        #plt.ylabel(r"$\tau$")
        
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.ylim([1e4,1e5])
        #plt.tight_layout()

        plt.show()

    if name == 'weissman_2010_2':
        plt.figure()
        plotted_times = [valley_crossing_times[N,sigma,mu,s,delta,rho] for rho in rho_range]
        plt.errorbar(rho_range, [np.average(x) for x in plotted_times], yerr=[np.std(x) for x in plotted_times], capsize=10)

        #theor = weissman_asexual_times(N, s, delta, mu_range)
        #plt.plot(rho_range, theor, ls = '--')
#        plt.xlabel('$\rho$')
#        plt.ylabel(r"$\tau$")
        
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()

        plt.show()
        
    if name == 'weissman_2010_3':
        plt.figure()
        plotted_times = [valley_crossing_times[N,sigma,mu,s,delta,rho] for rho in rho_range]
        plt.errorbar(rho_range, [np.average(x) for x in plotted_times], yerr=[np.std(x) for x in plotted_times], capsize=10)

        #theor = weissman_asexual_times(N, s, delta, mu_range)
        #plt.plot(rho_range, theor, ls = '--')
#        plt.xlabel('$\rho$')
#        plt.ylabel(r"$\tau$")
        
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()

        plt.show()

    if name == 'weissman_2010_4':
        plt.figure()
        plotted_times = [valley_crossing_times[N,sigma,mu,s,delta,rho] for rho in rho_range]
        plt.errorbar(rho_range, [np.average(x) for x in plotted_times], yerr=[np.std(x) for x in plotted_times], capsize=10)

        #theor = weissman_asexual_times(N, s, delta, mu_range)
        #plt.plot(rho_range, theor, ls = '--')
#        plt.xlabel('$\rho$')
#        plt.ylabel(r"$\tau$")
        
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()

        plt.show()


        
    for key in valley_crossing_times.keys():
        print key, len(valley_crossing_times[key])