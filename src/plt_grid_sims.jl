#!/usr/bin/env julia

## plt_grid_sims.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Make "grid" plots of crossing time versus δ and N.

using PyPlot, JLD, Glob, WaveFit, LaTeXStrings
#gr()

function load_and_plot(namestring::AbstractString)

    files = glob("output/vc_sims/grid_tunnel_*.jld")
    # glob all the output files

    # keys for this dict will be K value and s_u value

    for file in files
        f = load(file)
        println(f["params"])
    end
end
#sigmalist, UL_list, nt_results = load_and_plot("neut_tunnel")
#UL_vals = unique(UL_list)
#display(
load_and_plot("grid_tunnel")
