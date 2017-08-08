import numpy as np
import fwdpy11
import fwdpy11.model_params
import fwdpy11.fitness
import fwdpy11.wright_fisher

# plotting stuff
import matplotlib.pyplot as plt
import seaborn as sns
import powerlaw as pl

class BubbleRecorder(object):
    def __init__(self,ngens,N):
        #An index variable referring to where we will next store data
        self.i = 0
        self.weights = dict()
        
    def __call__(self,pop):
        for mut in range(len(pop.mutations)):
            # infinite sites model so each mutatio has unique position
            mkey = pop.mutations[mut].pos
            # if mutation neutral (neutral ones draft) 
            if pop.mutations[mut].neutral:
                # bubble is open
                if pop.mcounts[mut] > 0:
                    # update bubble size of currently segregating mutation
                    if mkey in self.weights:
                        self.weights[mkey] += pop.mcounts[mut]
                    # new segregating mutation; open bubble
                    else:
                        self.weights[mkey] = pop.mcounts[mut]
                    
def evolve_draft(ngens, N, s, mu, r, burn = 500, seed = 42):
    """
    Evolve drafting region and record bubble sizes
    """

    pop = fwdpy11.SlocusPop(N)
    rng=fwdpy11.GSLrng(seed)
    
    params = fwdpy11.model_params.SlocusParams(
        nregions = [fwdpy11.Region(0, 1, 1)],
        sregions = [fwdpy11.ExpS(0, 1, 1, s, 1.0)],
        recregions = [fwdpy11.Region(0, 1, 1)],
        gvalue = fwdpy11.fitness.SlocusAdditive(2.0),
        demography = np.array([N]*ngens, dtype=np.uint32),
        rates=(mu, mu, r)
        )
    recorder = BubbleRecorder(ngens, N)
    
    fwdpy11.wright_fisher.evolve(rng,pop,params,recorder)

    return pop, np.array(list(recorder.weights.values()))


# plot using the 'powerlaw' package and show fit
def plot_weight_powerlaw(weights, fit = None, loc = 3):
    if fit == None:
        fit = pl.Fit(weights, discrete=True)

    fit.plot_pdf(color = 'b', label = 'empirical')
    fit.power_law.plot_pdf(color = 'b', linestyle = '--',
                               label=r'power law fit: $\alpha$ = ' + str(np.round(fit.alpha, 3)))
    plt.legend(loc=3)

    return fit

# plot 'by hand'. borrow from TK code
def plot_weight_dist(weights, kde = True, hist = False, bins = 20, loc = 3):

    lmin = np.log10(min(weights))
    lmax = np.log10(max(weights))
    logbins = np.logspace(lmin, lmax, bins)
    h = np.histogram(weights, bins = logbins)
    
    w_bins = np.logspace(lmin, lmax, bins)
    w_2 = w_bins**-2
    w_32 = w_bins**-1.5
    w_1 = w_bins**-1

    plt.plot(h[1][:-1], h[0], label = 'empirical')
    plt.plot(w_bins, w_2, label = '$w^{-2}$')
    plt.plot(w_bins, w_32, label = '$w^{-3/2}$')
    plt.plot(w_bins, w_1, label = '$w^{-1}$')

    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc = loc)
    plt.show()

### ipython code: this produces a very good fit to the drift case of -alpha = -3/2
# %time pop, weights = evolve_draft(20000, 10000, 0.01, 0.001, 0.001, seed=314)
# plot_weight_powerlaw(weights)
