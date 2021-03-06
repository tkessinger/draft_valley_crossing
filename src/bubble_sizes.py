#!/usr/bin/env python
'''
simulate_population.py
Author: Taylor Kessinger
Date: June 10, 2017
Description: Generate bubble size distribution P(w) for a single mutant using FFPopSim.
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
import os
import cPickle as pickle

import seaborn as sns


t0 = time()

parser = argparse.ArgumentParser()
parser.add_argument('--N', type=int, default=1500)
parser.add_argument('--s', type=float, default=0.1)
parser.add_argument('--delta', type=float, default=0.0)
parser.add_argument('--rho', type=float, default=0.0001)
parser.add_argument('--s_mean', type=float, default=0)
parser.add_argument('--sigma', type=float, default=0.05)
parser.add_argument('--mu', type=float, default=0.00001)
parser.add_argument('--numsamples',type=int,default=10)
parser.add_argument('--runno', type=int, default=0)
parser.add_argument('--runtime', type=int, default=10000)
parser.add_argument('--name', type=str, default='test')

args=parser.parse_args()

L = 1000 #number of loci. 1000 is about as many as we're likely to need.
N = args.N #population size. recall that Nsigma is the main parameter.
mu = args.mu #per site mutation rate for focal mutations
sigma = args.sigma #fitness variance
rho = args.rho #per generation outcrossing rate
cx = 0 #per locus crossover rate
#unfortunately, nonzero values of rho and cx increase the runtime somewhat

s = args.s
delta = args.delta

nothing = 1e-12

mut_loci = 10

mutant_loc = np.array([np.int(L*(x+0.5)/mut_loci) for x in range(mut_loci)])
#mutant_loc = np.int(L*0.5) #locations of the focal mutations: this should avoid "edge" effects
#i am not sure why int32 is required

deleterious_effects = [-delta/2]
#surprisingly, these are the coefficients needed for single and double mutants to have the desired fitnesses.
#FFPopSim works on a {-1,1} basis, not a {0,1} basis, so
#F = \sum_i f_i s_i + \sum_{i<j} f_{ij} s_i s_j + ...
#the multiplication of indices is what makes the above necessary

pop = h.haploid_highd(L, number_of_traits=2, number_of_mutant_loci = 1, most_polymorphic = True) #high-dimensional simulation
#note: h.haploid_lowd() implements epistasis more easily, but we need highd to study draft
#note: all_polymorphic is an infinite sites model.
#loci whose frequencies hit 0 or 1 are automatically "injected" somewhere.
#if they fixed, their selection coefficients are flipped.
#most_polymorphic is just a hacked-together version that allows us to control the mutation rates of some loci.
pop.carrying_capacity = N
pop.mutation_rate = nothing
#pop._set_mutant_loci([mutant_loc])
pop._set_mutant_loci(mutant_loc)

pop.recombination_model = h.FREE_RECOMBINATION #recombination model, can be either FREE_RECOMBINATION or CROSSOVERS
pop.outcrossing_rate = rho
pop.crossover_rate = cx

pop.set_wildtype(N)

if args.s_mean:
    s_mean = args.s_mean #average effect of beneficial mutations: this turns out to be mostly irrelevant due to sigma, above
else:
    s_mean = 1e-8
base_burn_time = 0.1*N #number of generations to wait for population to equilibrate: ideally this should be of order T_MRCA
burn_time = base_burn_time
burned_in = False #flag for whether the population has equilibrated
run_time = args.runtime
#run_time = np.max([1e4,10*N]) #number of generations to simulate

num_attempts = 0 #attempts at crossing the valley, i.e., number of deleterious or neutral single mutants
mutation_in_play = False #flag to see whether the mutations are currently segregating
weights = []
bubble_openings = []
bubble_closings = []

dm_successful_crossings = []

mut_freqs = []


selection_coefficients = np.random.exponential(s_mean, size = L) #initialize fitness function
#selection_coefficients = np.array([s_mean for x in range(L)])
selection_coefficients[mutant_loc]=1e-10
#NOTE: due to an FFPopSim quirk, selection coefficients of 0 are incompatible with all_polymorphic.
#use something like 1e-10 instead.

pop.set_trait_additive(selection_coefficients,0)
pop.set_trait_additive(np.zeros(L),1)

#pop.trait_weights[0] = 0
#print pop.trait_weights

for loc in mutant_loc:
    pop.add_trait_coefficient(deleterious_effects[0],loc,1)
#mutant_loc might need to be a list
afreqs = []

if not args.s_mean:
    prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_delta_'+str(delta)+'_sigma_'+str(sigma)
else:
    prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_delta_'+str(delta)+'_sigma_'+str(s_mean)
if not os.path.exists('output/bubble_size_sims_'+args.name):
    os.makedirs('output/bubble_size_sims_'+args.name)

vars = []

while pop.generation < run_time:
#    if pop.get_allele_frequency(mutant_loc) == 0:
#        pop.flip_single_locus(mutant_loc)
    #print pop.generation
    pfit = pop.get_trait_statistics(0)
    #pfit1 = pop.get_trait_statistics(1)
    #fitvar = pop.get_fitness_statistics().variance
    #vars.append(pfit.variance)
    #print np.sqrt(fitvar), np.sqrt(pfit.variance), np.sqrt(pfit1.variance), sigma, pop.trait_weights
    #pfit = pop.get_fitness_statistics()
    if pfit.variance > 0 and pop.generation > np.min([100, burn_time]) and not args.s_mean:
        multiplier = sigma/np.sqrt(pfit.variance)
        #pop.trait_weights -= [pop.trait_weights[0],0]
        #pop.trait_weights += [multiplier,0]
        #print pop.get_trait_additive(0)[0]
        #pop.trait_weights[1] = multiplier
        #print multiplier, pop.trait_weights
        selection_coefficients *= multiplier
        selection_coefficients[mutant_loc]=1e-10
        pop.set_trait_additive(selection_coefficients,0)
        #normalizes the selection coefficients so that sigma stays constant.
        #the above is clunky compared to pop.trait_weights *= np.array([multiplier,1]),
        #but FFPopSim seems not to like trait weights higher than 1,
        #and it will gladly let trait weights collapse to 0 when there is more than one trait
    if pfit.variance == 0 and pop.generation > 100:
        #print "error: zero variance"
        pass
    if pop.generation > burn_time and burned_in == False: #if the population has equilibrated
        pop.trait_weights += np.array([0.0,1.0]) #"turn on" the focal loci
        #print "traits updated", pop.trait_weights
        pop.mutation_rate = mu #set the focal loci to start mutating
        burned_in = True
    if burned_in == True:
        if pop.get_allele_frequency(mutant_loc) != 0 and mutation_in_play == False:
            #if either mutant has appeared
            mutation_in_play = True
            bubble_openings.append(pop.generation)
            num_attempts += 1
            #print pop.get_allele_frequency(mutant_loc)*pop.N, int(np.ceil(pop.get_allele_frequency(mutant_loc)*pop.N))
            weights.append([int(np.ceil(pop.get_allele_frequency(mutant_loc)*pop.N))]) #start a new bubble of size 1
            pop.mutation_rate = nothing #turn off the bubbles
        elif mutation_in_play == True:
            #print pop.get_derived_allele_frequency(mutant_loc), pop.N*pop.get_allele_frequency(mutant_loc)
            if pop.get_allele_frequency(mutant_loc) == 0:
                mutation_in_play = False #if it's gone extinct, close the bubble
                #print "bubble closed"
                bubble_closings.append(pop.generation)
                pop.mutation_rate = mu
            elif pop.get_allele_frequency(mutant_loc) < 0.95:
                weights[-1].append(int(np.ceil(pop.get_allele_frequency(mutant_loc)*pop.N))) #add the current size to the bubble history
                    #ceil is needed to avoid rounding errors
            else:# pop.get_allele_frequency(mutant_loc) > 0.99*N:
                allele_freqs = pop.get_allele_frequencies() #if the locus fixes due to drift/draft
                allele_freqs[mutant_loc] = nothing #set it back to zero
                pop.set_allele_frequencies(allele_freqs,pop.N)
                weights[-1].append('fixed')
                mutation_in_play = False
                pop.mutation_rate = mu
                #print "removed rogue locus"
        #if not pop.generation%(10**4):
        '''
        if num_attempts > 0 and not (num_attempts%1e2) and mutation_in_play==False:
            #print "shuffling selection coefficients"
            #the below "clears out" the loci and reinitializes the population with the same appropriate allele frequencies
            allele_freqs = pop.get_allele_frequencies()
            allele_freqs[mutant_loc] = 0
            pop.set_allele_frequencies(allele_freqs,pop.N)
            pop.carrying_capacity = N #reset the carrying capacity to its desired value
            
            #selection_coefficients = np.random.exponential(s_mean, size = L) #we might as well change up the background fitnesses at this point
            selection_coefficients = np.array([s_mean for x in range(L)])

            selection_coefficients[mutant_loc]=1e-10
            
            if mutation_in_play:
                mutation_in_play = False
                weights[-1].append('bubble closed')
                #print "force closed bubble"
            
            num_attempts += 1
            pop.mutation_rate = nothing
            burn_time = pop.generation + base_burn_time #we might as well let the population equilibrate again for a bit
            #print burn_time
            burned_in = False
            pop.trait_weights += np.array([0,-1.0]) #"turn off" the focal loci: just a failsafe
        '''            

                
    mut_freqs.append(pop.get_allele_frequency(mutant_loc))
    afreqs.append([pop.get_derived_allele_frequency(L/4),pop.get_derived_allele_frequency(3*L/4)])
    if not pop.generation%(10**4):
        with open('output/bubble_size_sims_'+args.name+'/'+prefix+'_weights_'+str(args.runno)+'.pkl', 'w') as weights_file:
            pickle.dump(weights, weights_file)
        with open('output/bubble_size_sims_'+args.name+'/'+prefix+'_openings_'+str(args.runno)+'.pkl', 'w') as openings_file:
            pickle.dump(bubble_openings, openings_file)
        with open('output/bubble_size_sims_'+args.name+'/'+prefix+'_closings_'+str(args.runno)+'.pkl', 'w') as closings_file:
            pickle.dump(bubble_closings, closings_file)
    pop.evolve(1)

with open('output/bubble_size_sims_'+args.name+'/'+prefix+'_weights_'+str(args.runno)+'.pkl', 'w') as weights_file:
    pickle.dump(weights, weights_file)
with open('output/bubble_size_sims_'+args.name+'/'+prefix+'_openings_'+str(args.runno)+'.pkl', 'w') as openings_file:
    pickle.dump(bubble_openings, openings_file)
with open('output/bubble_size_sims_'+args.name+'/'+prefix+'_closings_'+str(args.runno)+'.pkl', 'w') as closings_file:
    pickle.dump(bubble_closings, closings_file)


if __name__ == '__main__':
    plt.figure()
    plt.plot(mut_freqs[::1])
    ax=plt.gca()
    ax.set_yscale('log')
    plt.show()
    afreqs = np.array(afreqs)
    plt.figure()
    plt.plot(afreqs[::1,0])
    plt.plot(afreqs[::1,1])
    plt.show()

    
    pruned_weights = []
    pruned_tots = []
    
    for weight in weights:
        if type(weight[-1]) != str:
            pruned_weights.append(weight)
            pruned_tots.append(weight)
    bubble_size_sums = [np.sum(weight) for weight in pruned_weights]
    bubble_size_tots = []
    [[bubble_size_tots.append(x) for x in weight] for weight in pruned_weights]
    plt.figure()
    w_theor_bins = np.logspace(2,7,50)
    w_2 = w_theor_bins**-2
    w_32 = w_theor_bins**-1.5
    w_1 = w_theor_bins**-1


    logbins = np.logspace(2,7,20)
    binwidths = logbins[1:] - logbins[:-1]
    #for sigma in [sigma_range[name][0],sigma_range[name][-1]]:
    #    for delta in [delta_range[name][0],delta_range[name][1]]:
    data = np.histogram(bubble_size_sums,bins=logbins,range=(10,1e4))#, density=True)
    hist_dat = 1.0*data[0]
    hist_dat*=w_32[0]/hist_dat[0]
    hist_dat /= (data[1][1:] - data[1][:-1])
    plt.plot(data[1][:-1],hist_dat,label='$\sigma$ = '+str(sigma) + ', $\delta$ = '+str(delta))
    
    plt.plot(w_theor_bins,w_2,label='$w^{-2}$')
    plt.plot(w_theor_bins,w_32,label='$w^{-3/2}$')
    plt.plot(w_theor_bins,w_1,label='$w^{-1}$')
    plt.legend(loc=1)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()
