#!/usr/bin/env python
'''
simulate_population.py
Author: Taylor Kessinger
Date: September 24, 2016
Description: Simulation of valley crossing using FFPopSim.
    Also TODO: direct calculation of Nsigma_b.
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


t0 = time()

parser = argparse.ArgumentParser()
parser.add_argument('--N', type=int, default=1000)
parser.add_argument('--s', type=float, default=0.1)
parser.add_argument('--delta', type=float, default=0.0)
parser.add_argument('--rho', type=float, default=0.0)
parser.add_argument('--sigma', type=float, default=0.01)
parser.add_argument('--mu', type=float, default=0.0001)
parser.add_argument('--numsamples',type=int,default=10)
parser.add_argument('--runno', type=int, default=0)
parser.add_argument('--runtime', type=int, default=100000)
parser.add_argument('--name', type=str, default='test')
parser.add_argument('--fixed_s', type=bool, default=False)

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

mutant_locs = np.array([np.int32(L*0.35), np.int32(L*0.65)]) #locations of the focal mutations: this should avoid "edge" effects
#i am not sure why int32 is required

deleterious_effects = [s/4,s/4]
beneficial_effect = delta/2+s/4
#surprisingly, these are the coefficients needed for single and double mutants to have the desired fitnesses.
#FFPopSim works on a {-1,1} basis, not a {0,1} basis, so
#F = \sum_i f_i s_i + \sum_{i<j} f_{ij} s_i s_j + ...
#the multiplication of indices is what makes the above necessary

pop = h.haploid_highd(L, number_of_traits=2, number_of_mutant_loci=2, most_polymorphic = True) #high-dimensional simulation
#note: h.haploid_lowd() implements epistasis more easily, but we need highd to study draft
#note: all_polymorphic is an infinite sites model.
#loci whose frequencies hit 0 or 1 are automatically "injected" somewhere.
#if they fixed, their selection coefficients are flipped.
#most_polymorphic is just a hacked-together version that allows us to control the mutation rates of some loci.
pop.carrying_capacity = N
pop.mutation_rate = 0
pop._set_mutant_loci(mutant_locs)
pop.recombination_model = h.FREE_RECOMBINATION #recombination model, can be either FREE_RECOMBINATION or CROSSOVERS
pop.outcrossing_rate = rho
pop.crossover_rate = cx

pop.set_wildtype(N)

if args.fixed_s:
    s_mean = sigma
    sigma_hist = []
else:
    s_mean = 1e-2 #average effect of beneficial mutations: this turns out to be mostly irrelevant due to sigma, above
base_burn_time = 0.1*N #number of generations to wait for population to equilibrate: ideally this should be of order T_MRCA
burn_time = base_burn_time
burned_in = False #flag for whether the population has equilibrated
#run_time = np.max([1e6,10*N]) #number of generations to simulate
run_time = args.runtime

num_attempts = 0 #attempts at crossing the valley, i.e., number of deleterious or neutral single mutants
num_beneficial = 0 #number of beneficial double mutants
num_successes = 0 #successful crossings
mutations_in_play = [False,False] #flag to see whether the mutations are currently segregating
double_mut_in_play = False
weights = [[],[]] #a record of the size of mutant bubbles
dm_weights = []
bubble_openings = [[],[]]
dm_bubble_openings = []
bubble_closings = [[],[]]
dm_bubble_closings = []

dm_successful_crossings = []

mut_freqs = [[],[]]
dm_freqs = []



selection_coefficients = np.random.exponential(s_mean, size = L) #initialize fitness function
selection_coefficients[mutant_locs]=1e-10
#NOTE: due to an FFPopSim quirk, selection coefficients of 0 are incompatible with all_polymorphic.
#use something like 1e-10 instead.

pop.set_trait_additive(selection_coefficients,0)
pop.set_trait_additive(np.zeros(L),1)

pop.add_trait_coefficient(deleterious_effects[0],[mutant_locs[0]],1)
pop.add_trait_coefficient(deleterious_effects[1],[mutant_locs[1]],1) #let the deleterious loci start drifting
pop.add_trait_coefficient(beneficial_effect, mutant_locs,1) #set the epistatic effect

afreqs = []

if args.fixed_s:
    prefix = 'N_'+str(N)+'_smean_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
else:
    prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)

if not os.path.exists('output/valley_crossing_sims_'+args.name):
    os.makedirs('output/valley_crossing_sims_'+args.name)

total_number_of_mutations = 0


while pop.generation < run_time:
    pfit = pop.get_trait_statistics(0)
    if pfit.variance > 0 and pop.generation > np.min([100, base_burn_time]) and not args.fixed_s:
        multiplier = sigma/np.sqrt(pfit.variance)
        pop.trait_weights -= [pop.trait_weights[0],0]
        pop.trait_weights += [multiplier,0]
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
        total_number_of_mutations += pop.number_of_mutations[-1]
        for di, del_locus in enumerate(mutant_locs):
            if pop.get_allele_frequency(del_locus) != 0 and mutations_in_play[di] == False:
            #if either mutant has appeared
                mutations_in_play[di] = True
                bubble_openings[di].append(pop.generation)
                num_attempts += 1
                weights[di].append([1]) #start a new bubble of size 1
            elif mutations_in_play[di] == True:
                if pop.get_allele_frequency(del_locus) == 0:
                    mutations_in_play[di] = False #if it's gone extinct, close the bubble
                    bubble_closings[di].append(pop.generation)
                else:
                    weights[di][-1].append(int(np.ceil(pop.get_allele_frequency(del_locus)*pop.N))) #add the current size to the bubble history
                    #ceil is needed to avoid rounding errors
        if pop.get_pair_frequency(mutant_locs[0], mutant_locs[1]) != 0 and double_mut_in_play == False:
        #if the double mutant has appeared
            double_mut_in_play = True
            dm_bubble_openings.append(pop.generation)
            num_beneficial += 1
            dm_weights.append([int(np.ceil(pop.get_pair_frequency(mutant_locs[0],mutant_locs[1])*pop.N))])
            #start a double mutant bubble of appropriate size (again, ceil avoids rounding errors)
        elif double_mut_in_play == True:
            if pop.get_pair_frequency(mutant_locs[0],mutant_locs[1]) == 0:
                dm_bubble_closings.append(pop.generation)
                double_mut_in_play = False
            elif pop.get_pair_frequency(mutant_locs[0],mutant_locs[1]) > 0.5: #this is a decent proxy for establishment
                #frequency 0.9 seems like about as good a cutoff as any
                dm_successful_crossings.append(pop.generation-burn_time) #record the amount of time elapsed since the last burn-in
                print 'successful crossing', pop.generation
                
                #the below "clears out" the loci and reinitializes the population with the same appropriate allele frequencies
                allele_freqs = pop.get_allele_frequencies()
                allele_freqs[mutant_locs] = 0
                pop.set_allele_frequencies(allele_freqs,pop.N)
                pop.carrying_capacity = N #reset the carrying capacity to its desired value
                
                #selection_coefficients = np.random.exponential(s_mean, size = L) #we might as well change up the background fitnesses at this point
                #selection_coefficients[mutant_locs]=1e-10
                
                double_mut_in_play = False
                
                pop.mutation_rate = 0
                burn_time = pop.generation + base_burn_time #we might as well let the population equilibrate again for a bit
                burned_in = False
                pop.trait_weights += np.array([0,-1.0]) #"turn off" the focal loci: just a failsafe

                
            else:
                dm_weights[-1].append(int(np.ceil(pop.get_pair_frequency(mutant_locs[0],mutant_locs[1])*pop.N)))
    dm_freqs.append(pop.get_pair_frequency(mutant_locs[0],mutant_locs[1]))
    mut_freqs[0].append(pop.get_allele_frequency(mutant_locs[0]))
    mut_freqs[1].append(pop.get_allele_frequency(mutant_locs[1]))
    if not pop.generation%(10**2) and pop.generation > base_burn_time and args.fixed_s:
        sigma_hist.append([np.sqrt(pop.get_trait_statistics(0).variance), np.sqrt(pop.get_fitness_statistics().variance)])
    if not pop.generation%(10**4):
        print pop.generation
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_dm_weights_'+str(args.runno)+'.pkl', 'w') as dm_weights_file:
            pickle.dump(dm_weights, dm_weights_file)
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_weights_'+str(args.runno)+'.pkl', 'w') as weights_file:
            pickle.dump(weights, weights_file)
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_dm_openings_'+str(args.runno)+'.pkl', 'w') as dm_openings_file:
            pickle.dump(bubble_openings, dm_openings_file)
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_openings_'+str(args.runno)+'.pkl', 'w') as openings_file:
            pickle.dump(dm_bubble_openings, openings_file)
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_dm_closings_'+str(args.runno)+'.pkl', 'w') as dm_closings_file:
            pickle.dump(dm_bubble_closings, dm_closings_file)
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_closings_'+str(args.runno)+'.pkl', 'w') as closings_file:
            pickle.dump(bubble_closings, closings_file)        
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_dm_successes_'+str(args.runno)+'.pkl', 'w') as dm_successes_file:
            pickle.dump(dm_successful_crossings, dm_successes_file)
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_effective_mut_rate_'+str(args.runno)+'.pkl', 'w') as mut_rate_file:
            pickle.dump([total_number_of_mutations, pop.generation - min([100,base_burn_time])], mut_rate_file)
        if args.fixed_s:
            with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_sigma_hist_'+str(args.runno)+'.pkl', 'w') as sigma_hist_file:
                pickle.dump(sigma_hist, sigma_hist_file)
        with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_current_generation.dat', 'w') as gen_file:
            gen_file.write(str(pop.generation))
    pop.evolve(1)

with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_dm_weights_'+str(args.runno)+'.pkl', 'w') as dm_weights_file:
    pickle.dump(dm_weights, dm_weights_file)
with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_weights_'+str(args.runno)+'.pkl', 'w') as weights_file:
    pickle.dump(weights, weights_file)
with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_dm_openings_'+str(args.runno)+'.pkl', 'w') as dm_openings_file:
    pickle.dump(bubble_openings, dm_openings_file)
with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_openings_'+str(args.runno)+'.pkl', 'w') as openings_file:
    pickle.dump(dm_bubble_openings, openings_file)
with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_dm_closings_'+str(args.runno)+'.pkl', 'w') as dm_closings_file:
    pickle.dump(dm_bubble_closings, dm_closings_file)
with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_closings_'+str(args.runno)+'.pkl', 'w') as closings_file:
    pickle.dump(bubble_closings, closings_file)
with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_dm_successes_'+str(args.runno)+'.pkl', 'w') as dm_successes_file:
    pickle.dump(dm_successful_crossings, dm_successes_file)
with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_effective_mut_rate_'+str(args.runno)+'.pkl', 'w') as mut_rate_file:
    pickle.dump([total_number_of_mutations, pop.generation - min([100,base_burn_time])], mut_rate_file)
if args.fixed_s:
    with open('output/valley_crossing_sims_'+args.name+'/'+prefix+'_sigma_hist_'+str(args.runno)+'.pkl', 'w') as sigma_hist_file:
        pickle.dump(sigma_hist, sigma_hist_file)


'''
if __name__ == '__main__':
    plt.figure()
    plt.plot(mut_freqs[0][::10])
    plt.plot(mut_freqs[1][::10])
    plt.plot(dm_freqs[::10])
    plt.show()
'''


#afreqs = np.array(afreqs)
#still needed:
#calculate Nsigma_b directly.

#print len(dm_weights)
#print len(weights[0]),len(weights[1])
