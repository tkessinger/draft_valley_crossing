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
parser.add_argument('--delta', type=float, default=0.0)
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

pop = h.haploid_highd(L, all_polymorphic=True) #high-dimensional simulation
#note: h.haploid_lowd() implements epistasis more easily, but we need highd to study draft
#note: all_polymorphic is an infinite sites model.
#loci whose frequencies hit 0 or 1 are automatically "injected" somewhere.
#if they fixed, their selection coefficients are flipped.
#most_polymorphic is just a hacked-together version that allows us to control the mutation rates of some loci.
pop.carrying_capacity = N
#pop.mutation_rate = nothing
#pop._set_mutant_loci([mutant_loc])
pop.recombination_model = h.FREE_RECOMBINATION #recombination model, can be either FREE_RECOMBINATION or CROSSOVERS
pop.outcrossing_rate = rho

pop.set_wildtype(N)

s_mean = 1e-6 #average effect of beneficial mutations: this turns out to be mostly irrelevant due to sigma, above
base_burn_time = 0.1*N #number of generations to wait for population to equilibrate: ideally this should be of order T_MRCA
burn_time = base_burn_time
burned_in = False #flag for whether the population has equilibrated
run_time = args.runtime
#run_time = np.max([1e6,10*N]) #number of generations to simulate


selection_coefficients = np.random.exponential(s_mean, size = L) #initialize fitness function
selection_coefficients[mutant_loc]=1e-10
#NOTE: due to an FFPopSim quirk, selection coefficients of 0 are incompatible with all_polymorphic.
#use something like 1e-10 instead.

pop.set_trait_additive(selection_coefficients,0)

prefix = 'N_'+str(N)+'_mu_'+str(mu)+'_s_'+str(s)+'_sigma_'+str(sigma)


while pop.generation < run_time:
#    if pop.get_allele_frequency(mutant_loc) == 0:
#        pop.flip_single_locus(mutant_loc)
    #print pop.generation
    pfit = pop.get_fitness_statistics()
    if pfit.variance > 0 and pop.generation > np.min([100, burn_time]):
        multiplier = sigma/np.sqrt(pfit.variance)
        selection_coefficients *= multiplier
        pop.set_trait_additive(selection_coefficients,0)
        #normalizes the selection coefficients so that sigma stays constant.
        #the above is clunky compared to pop.trait_weights *= np.array([multiplier,1]),
        #but FFPopSim seems not to like trait weights higher than 1,
        #and it will gladly let trait weights collapse to 0 when there is more than one trait
    if pfit.variance == 0 and pop.generation > 100:
        print "error: zero variance"
    if not pop.generation%(args.runtime/10):
        plt.figure()
        plt.hist(pop.get_traits()[:,0])
        ax=plt.gca()
        ax.set_yscale('log')
        plt.show()
    pop.evolve(1)
                



#afreqs = np.array(afreqs)
#still needed:
#calculate Nsigma_b directly.

#print len(dm_weights)
#print len(weights[0]),len(weights[1])
