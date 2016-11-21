#!/usr/bin/env python
'''
simulate_population.py
Author: Taylor Kessinger
Date: September 24, 2016
Description: Simulation of valley crossing using FFPopSim.
    Also TODO: direct calculation of Nsigma_b.
'''
import sys
import numpy as np
import FFPopSim as h
import random as rand
import matplotlib.pyplot as plt
from time import time
import argparse
import math
import cPickle as pickle

t0 = time()

parser = argparse.ArgumentParser()
parser.add_argument('--N', type=int, default=10000)
parser.add_argument('--s', type=float, default=0.1)
parser.add_argument('--delta', type=float, default=0.1)
parser.add_argument('--rho', type=float, default=0.0)
parser.add_argument('--sigma', type=float, default=0.1)
parser.add_argument('--mu', type=float, default=0.001)
parser.add_argument('--numsamples',type=int,default=10)
parser.add_argument('--runno', type=int, default=0)
parser.add_argument('--name', type=str, default='test')

args=parser.parse_args()

L = 300 #number of loci. 1000 is about as many as we're likely to need.
N = args.N #population size. recall that Nsigma is the main parameter.
mu = args.mu #per site mutation rate for focal mutations
sigma = args.sigma #fitness variance
rho = args.rho #per generation outcrossing rate
cx = 0 #per locus crossover rate
#unfortunately, nonzero values of rho and cx increase the runtime somewhat

s = args.s
delta = args.delta

fitness_check = N*10 #for debugging: number of generations to let pass before checking the gaussianity of the fitness distribution

mutant_locs = np.array([np.int32(L*0.35), np.int32(L*0.65)]) #locations of the focal mutations: this should avoid "edge" effects
#i am not sure why int32 is required
deleterious_effects = [delta,delta] #fitness costs of single mutants
beneficial_effect = s #fitness advantage of double mutant

pop = h.haploid_highd(L, most_polymorphic = True) #high-dimensional simulation
#note: h.haploid_lowd() implements epistasis more easily, but we need highd to study draft
#note: all_polymorphic is an infinite sites model.
#loci whose frequencies hit 0 or 1 are automatically "injected" somewhere.
#if they fixed, their selection coefficients are flipped.
#i believe FFPopSim has some kind of flag to see whether fixation/extinction have occurred (need to look for this or write it).
pop.carrying_capacity = N
pop.mutation_rate = mu
pop._set_mutant_loci(mutant_locs)
pop.recombination_model = h.CROSSOVERS #recombination model, can be either FREE_RECOMBINATION or CROSSOVERS
pop.outcrossing_rate = rho
pop.crossover_rate = cx

pop.set_wildtype(N)

s_mean = 0.01 #average effect of beneficial mutations: this turns out to be mostly irrelevant due to sigma, above
burn_time = N #number of generations to wait for population to equilibrate: ideally this should be of order T_MRCA
burned_in = False #flag for whether the population has equilibrated
run_time = 5*N #number of generations to simulate

num_attempts = 0 #attempts at crossing the valley, i.e., number of deleterious or neutral single mutants
num_beneficial = 0 #number of beneficial double mutants
num_successes = 0 #successful crossings
mutations_in_play = [False,False] #flag to see whether the mutations are currently segregating
double_mut_in_play = False
weights = [[],[]] #a record of the size of mutant bubbles
dm_weights = []

selection_coefficients = np.random.exponential(s_mean, size = L) #initialize fitness function
#NOTE: due to an FFPopSim quirk, selection coefficients of 0 are incompatible with all_polymorphic.
#use something like 1e-10 instead.

pop.set_trait_additive(selection_coefficients)
pop.add_trait_coefficient(-0.5, [mutant_locs[0]])
pop.add_trait_coefficient(-0.5, [mutant_locs[1]]) #set the fitness of the deleterious mutants to (effectively) zero.
pop.add_trait_coefficient(beneficial_effect, mutant_locs) #set the epistatic effect

afreqs = []

while pop.generation < run_time:
    if not pop.generation%fitness_check:
        #check the histogram for gaussianity manually
        pfhist = pop.get_fitness_histogram(100)
        plt.figure()
        plt.plot((pfhist[1][1:]+pfhist[1][:-1])*1.0/2,pfhist[0])
        ax=plt.gca()
        ax.set_yscale('log') #on a log scale, the fitness distribution should be parabolic
        plt.show()
    pfit = pop.get_fitness_statistics()
    if pfit.variance > 0 and pop.generation > 100:
        pop.trait_weights *= sigma/np.sqrt(pfit.variance) #normalize the selection coefficients so that sigma is constant.
        #TODO: exclude the focal loci from this calculation.
    if pop.generation > burn_time and burned_in == False: #if the population has equilibrated
        pop.add_trait_coefficient(deleterious_effects[0],[mutant_locs[0]])
        pop.add_trait_coefficient(deleterious_effects[1],[mutant_locs[1]]) #let the deleterious loci start drifting
        burned_in = True
    if burned_in == True:
        #this loop is currently never called due to a quirk in all_polymorphic (TODO).
        #in principle the idea is to see if the new locus has entered the population anew.
        #because all_polymorphic automatically seeds extinct alleles in one individual,
        #this is impossible to determine (frequency is never allowed to hit zero).
        #there are a few possible ways around this.
        for di, del_locus in enumerate(mutant_locs):
            if pop.get_allele_frequency(del_locus) > 0 and mutations_in_play[di] == False:
            #if either mutant has appeared
                mutations_in_play[di] = True
                num_attempts += 1
                weights[di].append([1]) #start a new bubble of size 1
            elif mutations_in_play[di] == True:
                if pop.get_allele_frequency(del_locus) == 0:
                    mutations_in_play[di] = False #if it's gone extinct, close the bubble
                else:
                    weights[di][-1].append(int(pop.get_allele_frequency(del_locus)*pop.N)) #add the current size to the bubble history
        if pop.get_pair_frequency(mutant_locs[0], mutant_locs[1]) > 0 and double_mut_in_play == False:
        #if the double mutant has appeared
            double_mut_in_play = True
            num_beneficial += 1
            dm_weights.append([1]) #start a double mutant bubble of size 1
        elif double_mut_in_play == True:
            if pop.get_pair_frequency(mutant_locs[0],mutant_locs[1]) == 0:
                double_mut_in_play = False
            else:
                dm_weights[-1].append(int(pop.get_pair_frequency(mutant_locs[0],mutant_locs[1])*pop.N))
    afreqs.append(pop.get_allele_frequencies())
    pop.evolve(1)


afreqs = np.array(afreqs)
#still needed:
#calculate Nsigma_b directly.
