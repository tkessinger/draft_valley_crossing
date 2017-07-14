#!/usr/bin/env python
'''
fixation_probability.py
Author: Taylor Kessinger
Date: June 10, 2017
Description: Generate fixation probabilities for a single mutant with fixed sigma.
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
import os
import cPickle as pickle


t0 = time()

parser = argparse.ArgumentParser()
parser.add_argument('--N', type=int, default=100)
parser.add_argument('--s', type=float, default=0.1)
parser.add_argument('--rho', type=float, default=0.0)
parser.add_argument('--sigma', type=float, default=0.01)
parser.add_argument('--mu', type=float, default=0.00001)
parser.add_argument('--numsamples',type=int,default=10)
parser.add_argument('--runno', type=int, default=0)
parser.add_argument('--name', type=str, default='test')
parser.add_argument('--runtime', type=int, default=100000)

args=parser.parse_args()

L = 200 #number of loci. 1000 is about as many as we're likely to need.
N = args.N #population size. recall that Nsigma is the main parameter.
mu = args.mu #per site mutation rate for focal mutations
sigma = args.sigma #fitness variance
rho = args.rho #per generation outcrossing rate
cx = 0 #per locus crossover rate
#unfortunately, nonzero values of rho and cx increase the runtime somewhat

s = args.s
delta = args.delta

nothing = 1e-12

mutant_loc = np.int(L*0.5) #locations of the focal mutations: this should avoid "edge" effects
#i am not sure why int32 is required

mutation_effect = [s/2]
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
pop._set_mutant_loci([mutant_loc])
pop.recombination_model = h.FREE_RECOMBINATION #recombination model, can be either FREE_RECOMBINATION or CROSSOVERS
pop.outcrossing_rate = rho
#pop.crossover_rate = cx

pop.set_wildtype(N)

s_mean = 1e-6 #average effect of beneficial mutations: this turns out to be mostly irrelevant due to sigma, above
base_burn_time = 0.1*N #number of generations to wait for population to equilibrate: ideally this should be of order T_MRCA
burn_time = base_burn_time
burned_in = False #flag for whether the population has equilibrated
run_time = args.runtime
#run_time = np.max([1e6,10*N]) #number of generations to simulate

num_attempts = 0 #attempts at crossing the valley, i.e., number of deleterious or neutral single mutants
mutation_in_play = False #flag to see whether the mutations are currently segregating
successful_crossings = 0
failed_crossings = 0

mut_freqs = []


selection_coefficients = np.random.exponential(s_mean, size = L) #initialize fitness function
selection_coefficients[mutant_loc]=1e-10
#NOTE: due to an FFPopSim quirk, selection coefficients of 0 are incompatible with all_polymorphic.
#use something like 1e-10 instead.

pop.set_trait_additive(selection_coefficients,0)
pop.set_trait_additive(np.zeros(L),1)

pop.add_trait_coefficient(mutation_effect[0],[mutant_loc],1)
#mutant_loc might need to be a list
afreqs = []

prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_s_'+str(s)+'_sigma_'+str(sigma)

if not os.path.exists('output/fix_prob_sims_'+args.name):
    os.makedirs('output/fix_prob_sims_'+args.name)


while pop.generation < run_time:
#    if pop.get_allele_frequency(mutant_loc) == 0:
#        pop.flip_single_locus(mutant_loc)
    #print pop.generation
    #pfit = pop.get_trait_statistics(0)
    pfit = pop.get_trait_statistics(0)
    if pfit.variance > 0 and pop.generation > 100:
        multiplier = sigma/np.sqrt(pfit.variance)
        selection_coefficients *= multiplier
        selection_coefficients[mutant_loc]=1e-10
        pop.set_trait_additive(selection_coefficients,0)
        #normalizes the selection coefficients so that sigma stays constant.
        #the above is clunky compared to pop.trait_weights *= np.array([multiplier,1]),
        #but FFPopSim seems not to like trait weights higher than 1,
        #and it will gladly let trait weights collapse to 0 when there is more than one trait
    if pfit.variance == 0 and pop.generation > 100:
        print "error: zero variance"
    if pop.generation > burn_time and burned_in == False: #if the population has equilibrated
        pop.trait_weights += np.array([0,1.0]) #"turn on" the focal loci
        pop.mutation_rate = mu #set the focal loci to start mutating
        burned_in = True
    if burned_in == True:
        if pop.get_allele_frequency(mutant_loc) != 0 and mutation_in_play == False:
            #if either mutant has appeared
            mutation_in_play = True
            num_attempts += 1
            pop.mutation_rate = nothing #turn off the bubbles
        elif mutation_in_play == True:
            #print pop.get_derived_allele_frequency(mutant_loc), pop.N*pop.get_allele_frequency(mutant_loc)
            if pop.get_allele_frequency(mutant_loc) == 0:
                mutation_in_play = False #if it's gone extinct, close the bubble
                #print "bubble closed"
                failed_crossings += 1
                pop.mutation_rate = mu
                    #ceil is needed to avoid rounding errors


            elif pop.get_allele_frequency(mutant_loc) == 1:
                allele_freqs = pop.get_allele_frequencies() #if the locus fixes due to drift/draft
                allele_freqs[mutant_loc] = 0 #set it back to zero
                pop.set_allele_frequencies(allele_freqs,pop.N)
                #print "successful crossing"
                successful_crossings += 1
                mutation_in_play = False
                pop.mutation_rate = mu
                #print "removed rogue locus"
        #if not pop.generation%(10**4):
        if num_attempts > 0 and not (num_attempts%1e2) and mutation_in_play==False:
            #the below "clears out" the loci and reinitializes the population with the same appropriate allele frequencies
            allele_freqs = pop.get_allele_frequencies()
            allele_freqs[mutant_loc] = 0
            pop.set_allele_frequencies(allele_freqs,pop.N)
            pop.carrying_capacity = N #reset the carrying capacity to its desired value
            
            selection_coefficients = np.random.exponential(s_mean, size = L) #we might as well change up the background fitnesses at this point
            selection_coefficients[mutant_loc]=1e-10
            
            if mutation_in_play:
                mutation_in_play = False
            
            num_attempts += 1
            pop.mutation_rate = nothing
            burn_time = pop.generation + base_burn_time #we might as well let the population equilibrate again for a bit
            #print burn_time
            burned_in = False
            pop.trait_weights += np.array([0,-1.0]) #"turn off" the focal loci: just a failsafe
            
#     if not pop.generation%(args.runtime/10):
#         plt.figure()
#         plt.hist(pop.get_traits()[:,0])
#         ax=plt.gca()
#         ax.set_yscale('log')
#         plt.show()
                
    mut_freqs.append(pop.get_allele_frequency(mutant_loc))
    if not pop.generation%(10**4):
        with open('output/fix_prob_sims_'+args.name+'/'+prefix+'_weights_'+str(args.runno)+'.pkl', 'w') as prob_file:
            pickle.dump([num_attempts, failed_crossings, successful_crossings], prob_file)
    pop.evolve(1)

with open('output/fix_prob_sims_'+args.name+'/'+prefix+'_weights_'+str(args.runno)+'.pkl', 'w') as prob_file:
    pickle.dump([num_attempts, failed_crossings, successful_crossings], prob_file)


if __name__ == '__main__':
    print failed_crossings, successful_crossings, 1.0*successful_crossings/(failed_crossings+successful_crossings)
    
    plt.figure()
    plt.plot(mut_freqs[::10])
    plt.show()



#afreqs = np.array(afreqs)
#still needed:
#calculate Nsigma_b directly.

#print len(dm_weights)
#print len(weights[0]),len(weights[1])
