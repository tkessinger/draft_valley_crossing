#!/usr/bin/env python
'''
plt_asex_crossings.py
Author: Taylor Kessinger
Date: July 24, 2018
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
import os

import seaborn as sns
sns.set_style('white')
sns.set_context('paper')
sns.set(font="CMU Sans Serif", font_scale=1.2,style='white')#, style='ticks')


numruns = {}
sigma_range = {}
N_range = {}
s_range = {}
delta_range = {}
mu_range = {}
rho_range = {}

run_time = {}

N_range = [np.int(x) for x in np.logspace(3,6,7)]
#    N_range[name] = [10**3,3*10**3, 10**4,10**5]#,3*10**4,10**5,3*10**5,10**6] #range of population size values
sigma_range = [5*10**-6,5*10**-2]
#delta_range = [0.0,10**-4,10**-2] #range of deleterious intermediate coefficients
delta_range = [0.0, 10**-4.5,10**-4,10**-3.5,10**-3,10**-2.5,10**-2,10**-1.5]
s_range_sweep = [10**-2] #range of selection coefficients for the adaptation
s_range_valley = [10**-1]
valley_name = 'compare_ratios_high_mu'
if valley_name == 'compare_ratios_high_mu':
    sweeps_name = 'compare_ratios_sweeps_high_mu'
    mu_range = [5*10**-4]
elif valley_name == 'compare_ratios_low_mu':
    mu_range = [5*10**-6]
    sweeps_name = 'compare_ratios_sweeps_low_mu'

elif valley_name == 'compare_ratios_very_high_mu':
    sweeps_name = 'compare_ratios_sweeps_very_high_mu'

    mu_range = [10**-4]
else:
    sweeps_name= 'compare_ratios_sweeps_2'
    mu_range = [5*10**-5]
numruns = 3
run_time = 1000000
rho = 0.0



def get_alpha_ratios(valley_name, sweeps_name, mu_range):
    valley_crossing_times = {}
    sweep_times = {}

    for N in N_range:
        for sigma in sigma_range:
            for mu in mu_range:
                for delta in delta_range:
                    for s in s_range_valley:
                        params= (N,sigma,mu,delta,s)
                        if not params in valley_crossing_times.keys():
                            valley_crossing_times[params] = []
                        for runno in range(numruns):
                            prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                            if os.path.isfile('output/valley_crossing_sims_'+valley_name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl'):
                                with open('output/valley_crossing_sims_'+valley_name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                    times = pickle.load(dm_successes_file)
                                    for time in times:
                                        valley_crossing_times[params].append(time)
                        print 'valley', params, len(valley_crossing_times[params]), np.average(valley_crossing_times[params]), np.std(valley_crossing_times[params])
                for s in s_range_sweep:
                    delta = 0.0
                    params= (N,sigma,mu,s)
                    if not params in sweep_times.keys():
                        sweep_times[params] = []
                    for runno in range(numruns):
                        name = sweeps_name
                        prefix = 'N_'+str(N)+'_sigma_'+str(sigma)+'_mu_'+str(mu)+'_rho_'+str(rho)+'_s_'+str(s)+'_delta_'+str(delta)
                        if os.path.isfile('output/single_sweep_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl'):
                            with open('output/single_sweep_sims_'+name+'/'+prefix+'_dm_successes_'+str(runno)+'.pkl', 'r') as dm_successes_file:
                                times = pickle.load(dm_successes_file)
                                for time in times:
                                    sweep_times[params].append(time)
                    print 'sweep', params, len(sweep_times[params]), np.average(sweep_times[params])
                
    sigma = sigma_range[0]
    #mu = 5*10**-4
    array1 = [[np.average(valley_crossing_times[N,sigma,mu,delta,s_range_valley[0]]) for delta in delta_range[1:]] for N in N_range]
    array1 = np.array(array1)
    #mu = 5*10**-5
    array2 = [[np.average(sweep_times[N,sigma,mu,s_range_sweep[0]]) for x in range(7)] for N in N_range]
    array2 = np.array(array2)

    sigma = sigma_range[1]
    array3 = [[np.average(valley_crossing_times[N,sigma,mu,delta,s_range_valley[0]]) for delta in delta_range[1:]] for N in N_range]
    array3 = np.array(array3)
    array4 = [[np.average(sweep_times[N,sigma,mu,s_range_sweep[0]]) for x in range(7)] for N in N_range]
    array4 = np.array(array4)
    
    return (array3/array4)/(array1/array2)

# 
# yticks = [0,3,3.5,4,4.5,5,5.5,6,6.5]
# xticks = [0,-4.5,-4,-3.5,-3,-2.5,-2,-1.5]
# 
# yticks = [r'$0$'] + [r'$10^{' + str(x) + '}$' for x in [3,3.5,4,4.5,5,5.5,6]]
# xticks = [r'$0$',r'$0$'] + [r'$10^{' + str(x) + '}$' for x in [-4.5,-4,-3.5,-3,-2.5,-2]]
# 
# 
# plt.figure()
# plt.imshow(np.log10(array1),cmap='Greens',vmin=2,vmax=4)
# plt.xlabel(r'$\delta$')
# plt.ylabel(r'$\log_{10}(N)$')
# plt.title(r'$\tau_{\mathrm{valley}}, \sigma = 5\times 10^{-6}$')
# ax=plt.gca()
# ax.set_xticklabels(xticks)
# ax.set_yticklabels(yticks)
# col = plt.colorbar()
# col.set_label(r'$\log_{10}\tau_{\mathrm{valley}}$')
# plt.tight_layout()
# 
# plt.savefig('tau_valley_-6.pdf',bbox_inches='tight')
# 
# plt.figure()
# plt.imshow(np.log10(array2),cmap='Greens',vmin=2,vmax=4)
# plt.xlabel(r'$\delta$')
# plt.ylabel(r'$\log_{10}(N)$')
# plt.title(r'$\tau_{\mathrm{sweep}}, \sigma = 5\times 10^{-6}$')
# ax=plt.gca()
# ax.set_xticklabels(xticks)
# ax.set_yticklabels(yticks)
# col = plt.colorbar()
# col.set_label(r'$\log_{10}\tau_{\mathrm{sweep}}$')
# plt.tight_layout()
# 
# plt.savefig('tau_sweep_-6.pdf',bbox_inches='tight')
# 
# 
# plt.figure()
# plt.imshow(np.log10(array1/array2),cmap='Greens')
# plt.xlabel(r'$\delta$')
# plt.ylabel(r'$\log_{10}(N)$')
# plt.title(r'$\alpha = \tau_{\mathrm{valley}}/\tau_{\mathrm{sweep}}, \sigma = 5\times 10^{-6}$')
# ax=plt.gca()
# ax.set_xticklabels(xticks)
# ax.set_yticklabels(yticks)
# col = plt.colorbar()
# col.set_label(r'$\log_{10}\alpha$')
# plt.tight_layout()
# 
# plt.savefig('alpha_-6.pdf',bbox_inches='tight')
# 
# sigma = sigma_range[1]
# array3 = [[np.average(valley_crossing_times[N,sigma,mu,delta,s_range_valley[0]]) for delta in delta_range[1:]] for N in N_range]
# array3 = np.array(array3)
# array4 = [[np.average(sweep_times[N,sigma,mu,s_range_sweep[0]]) for x in range(7)] for N in N_range]
# array4 = np.array(array4)
# 
# plt.figure()
# plt.imshow(np.log10(array3),cmap='Greens',vmin=2,vmax=4)
# plt.xlabel(r'$\delta$')
# plt.ylabel(r'$\log_{10}(N)$')
# plt.title(r'$\tau_{\mathrm{valley}}, \sigma = 5\times 10^{-2}$')
# ax=plt.gca()
# ax.set_xticklabels(xticks)
# ax.set_yticklabels(yticks)
# col = plt.colorbar()
# col.set_label(r'$\log_{10}\tau_{\mathrm{valley}}$')
# plt.tight_layout()
# 
# plt.savefig('tau_valley_-2.pdf',bbox_inches='tight')
# 
# plt.figure()
# plt.imshow(np.log10(array4),cmap='Greens',vmin=2,vmax=4)
# plt.xlabel(r'$\delta$')
# plt.ylabel(r'$\log_{10}(N)$')
# plt.title(r'$\tau_{\mathrm{sweep}}, \sigma = 5\times 10^{-2}$')
# ax=plt.gca()
# ax.set_xticklabels(xticks)
# ax.set_yticklabels(yticks)
# col = plt.colorbar()
# col.set_label(r'$\log_{10}\tau_{\mathrm{sweep}}$')
# plt.tight_layout()
# plt.savefig('tau_sweep_-2.pdf',bbox_inches='tight')
# 
# plt.figure()
# plt.imshow(np.log10(array3/array4),cmap='Greens')
# plt.xlabel(r'$\delta$')
# plt.ylabel(r'$\log_{10}(N)$')
# plt.title(r'$\alpha = \tau_{\mathrm{valley}}/\tau_{\mathrm{sweep}}, \sigma = 5\times 10^{-2}$')
# ax=plt.gca()
# ax.set_xticklabels(xticks)
# ax.set_yticklabels(yticks)
# col = plt.colorbar()
# col.set_label(r'$\log_{10}\alpha$')
# plt.tight_layout()
# 
# plt.savefig('alpha_-2.pdf',bbox_inches='tight')
# 
# 
# plt.show()
# 
# plt.figure()
# plt.imshow(np.log10((array3/array4)/(array1/array2)),cmap='Greens')
# plt.xlabel(r'$\delta$')
# plt.ylabel(r'$\log_{10}(N)$')
# plt.title(r'$\alpha_{\mathrm{\sigma = 5\times 10^{-2}}}/\alpha_{\mathrm{\sigma = 5\times 10^{-6}}}$')
# ax=plt.gca()
# ax.set_xticklabels(xticks)
# ax.set_yticklabels(yticks)
# col = plt.colorbar()
# col.set_label(r'$\log_{10}\alpha_{\mathrm{\sigma = 5\times 10^{-2}}}/\alpha_{\mathrm{\sigma = 5\times 10^{-6}}}$')
# plt.show()
# plt.tight_layout()
# 
# plt.savefig('alpha_ratio.pdf',bbox_inches='tight')
# 

yticks = [r'$0$'] + [r'$10^{' + str(x) + '}$' for x in [3,3.5,4,4.5,5,5.5,6]]
xticks = [r'$0$',r'$0$'] + [r'$10^{' + str(x) + '}$' for x in [-4.5,-4,-3.5,-3,-2.5,-2]]

from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib
 
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero
 
    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }
 
    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)
 
    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])
 
    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
 
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
 
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
 
    return newcmap
 
orig_cmap = matplotlib.cm.RdBu
shifted_cmap = shiftedColorMap(orig_cmap, midpoint=(0.6/1.9), name='shifted')
 
from mpl_toolkits.axes_grid1 import ImageGrid
''' 
fig = plt.figure(figsize=(10.75, 5))
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,2),
                 axes_pad=0.15,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )
axes = grid
im = axes[0].imshow(np.log10(array1/array2),cmap=shifted_cmap,vmin=-.6,vmax=1.3)
axes[0].set_xticklabels(xticks)
axes[0].set_yticklabels(yticks)
axes[0].set_xlabel(r'$\delta$')
axes[0].set_ylabel(r'$N$')
axes[0].set_title(r'$\sigma_1 = 5 \times 10^{-6}$')
im = axes[1].imshow(np.log10(array3/array4),cmap=shifted_cmap,vmin=-.6,vmax=1.3)
axes[1].set_xticklabels(xticks)
axes[1].set_yticklabels(yticks)
axes[1].set_xlabel(r'$\delta$')
axes[1].set_ylabel(r'$N$')
axes[1].set_title(r'$\sigma_2 = 5 \times 10^{-2}$')
cb = axes[1].cax.colorbar(im)
axes[1].cax.toggle_label(True)
cb.set_label_text(r'$\log_{10} (\alpha)$')
#plt.tight_layout()    # Works, but may still require rect parameter to keep colorbar labels visible
plt.show()
'''
alpha1 = get_alpha_ratios('compare_ratios_2', 'compare_ratios_sweeps_2', [5*10**-5])
alpha2 = get_alpha_ratios('compare_ratios_very_high_mu', 'compare_ratios_sweeps_very_high_mu', [10**-4])
alpha3 = get_alpha_ratios('compare_ratios_high_mu', 'compare_ratios_sweeps_high_mu', [5*10**-4])

begin = -0.8
end = 0.2

orig_cmap = matplotlib.cm.BrBG
shifted_cmap = shiftedColorMap(orig_cmap, midpoint=(-1.0*begin/(-begin+end)), name='shifted')


fig = plt.figure(figsize=(14.25, 6))
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,3),
                 axes_pad=0.15,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )
axes = grid
im = axes[0].imshow(np.log10(alpha1),cmap=shifted_cmap,vmin=begin,vmax=end)
axes[0].set_xticklabels(xticks)
axes[0].set_yticklabels(yticks)
axes[0].set_xlabel(r'$\delta$')
axes[0].set_ylabel(r'$N$')
axes[0].set_title(r'$\mu_1 = 5 \times 10^{-5}$')
im = axes[1].imshow(np.log10(alpha2),cmap=shifted_cmap,vmin=begin,vmax=end)
axes[1].set_xticklabels(xticks)
axes[1].set_yticklabels(yticks)
axes[1].set_xlabel(r'$\delta$')
axes[1].set_ylabel(r'$N$')
axes[1].set_title(r'$\mu_2 = 10^{-4}$')
im = axes[2].imshow(np.log10(alpha3),cmap=shifted_cmap,vmin=begin,vmax=end)
axes[2].set_xticklabels(xticks)
axes[2].set_yticklabels(yticks)
axes[2].set_xlabel(r'$\delta$')
axes[2].set_ylabel(r'$N$')
axes[2].set_title(r'$\mu_3 = 5 \times 10^{-4}$')
#axes[0].set_title(r'$\sigma_1 = 5 \times 10^{-6}$')
cb = axes[0].cax.colorbar(im)
axes[0].cax.toggle_label(True)
cb.set_label_text(r'$\log_{10}(\alpha_{\sigma_2}/\alpha_{\sigma_1})$')
#plt.tight_layout()    # Works, but may still require rect paramater to keep colorbar labels visible
plt.show()


