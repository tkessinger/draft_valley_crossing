#!/usr/bin/env python3

import numpy as np
import fwdpy11
import fwdpy11.model_params
import fwdpy11.fitness
import fwdpy11.wright_fisher
import powerlaw as pl
import pickle
import fcntl
from argparse import ArgumentParser
from datetime import datetime

# pybind11 C++ recorder
import cppimport
cppimport.set_quiet(False)
bubble_recorder = cppimport.imp("bubble_recorder")


def evolve_draft(ngens, N, s, mu, r, burn=500, reps=10, seed=42):
    """
    Evolve drafting region and record bubble sizes
    """

    rng = fwdpy11.GSLrng(seed)

    pdict = dict(nregions=[fwdpy11.Region(0, 1, 1)],
                 sregions=[fwdpy11.ExpS(0, 1, 1, s, 1.0)],
                 recregions=[fwdpy11.Region(0, 1, 1)],
                 rates=(mu, mu, r))

    # shutup the warnings from the powerlaw fit
    np.seterr(divide='ignore', invalid='ignore')

    pops = list()
    weights = list()
    alpha = list()
    w_mean = list()
    w_var = list()
    for i in range(reps):
        pop = fwdpy11.SlocusPop(N)

        # burnin runs
        pdict.update(demography=np.array([N]*burn, dtype=np.uint32))
        params = fwdpy11.model_params.SlocusParams(**pdict)
        fwdpy11.wright_fisher.evolve(rng, pop, params)

        # recorded runs
        pdict.update(demography=np.array([N]*ngens, dtype=np.uint32))
        recorder = bubble_recorder.BubbleRecorder()
        fwdpy11.wright_fisher.evolve(rng, pop, params, recorder)

        wvec = np.array(list(recorder.weights.values()))
        pops.append(pop)
        weights.append(wvec)
        alpha.append(pl.Fit(wvec, discrete=True, verbose=False).alpha)
        w_mean.append(np.array(recorder.w_mean))
        w_var.append(np.array(recorder.w_var))

    return pops, weights, alpha, w_mean, w_var


# p, w, a, m, v = evolve_draft(300, 1000, 0.01, 0.001, 0.0001,
#                              seed=314, burn=100, reps=10)


def main():
    parser = ArgumentParser(description='mutant bubble weights')
    parser.add_argument("--N", type=int, default=10000, help="population size")
    parser.add_argument("--s", type=float, default=0.001, help="selection coefficient")
    parser.add_argument("--mu", type=float, default=0.0001, help="mutation rate")
    parser.add_argument("--r", type=float, default=0.00001, help="recombination rate")
    parser.add_argument("--ngens", type=int, default=10000, help="number of generations (after burnin)")
    parser.add_argument("--burn", type=int, default=5000, help="burnin generations")
    parser.add_argument("--reps", type=int, default=10, help="number of replicates")
    parser.add_argument("--seed", type=int, default=42, help="random number seed")
    parser.add_argument("--file", default='data.pkl', help="output file")

    args = parser.parse_args()

    start = datetime.now()
    print (':: sim started at ' + str(start))
    pops, weights, alpha, means, vars = evolve_draft(args.ngens,
                                                     args.N,
                                                     args.s,
                                                     args.mu,
                                                     args.r,
                                                     args.burn,
                                                     args.reps,
                                                     args.seed)

    with open(args.file, 'ab') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        pickle.dump(
            {'pars': dict(args._get_kwargs()),
             'pops': pops,
             'weights': weights,
             'alpha': alpha,
             'means': means,
             'vars': vars}, f)

    end = datetime.now
    print(':: sim finished at ' + str(end) + ' and took ' + str(end-start))


if __name__ == "__main__":
    main()
