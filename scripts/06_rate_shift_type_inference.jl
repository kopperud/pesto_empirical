## Inference of the simulated datasets

using Revise
using Pesto
using Distributions
using CSV
using JLD2
using ProgressMeter

#= 
These are the outputs I want for the inference:
    * Newick string
    * N matrix (3-dim)
    * Tip rates ??
    * Branch rates CSV file
    * Constant-rate Birth-Death Process estimates
    * Estimate for η
=#

## iter 208 was not a good look -- segfault or something

n_iters = 1000
prog = ProgressMeter.Progress(n_iters; desc = "Inference: ");
for i in 532:n_iters
    phy = readtree(string("data/simulations/rate_shift_type/", i ,".tre"))
    ρ = 1.0
    data = SSEdata(phy, ρ);

    λml, μml = estimate_constant_bdp(data)

    H = 0.587
    n = 5

    dλ = LogNormal(log(λml), H)
    dμ = LogNormal(log(µml), H)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)

    ηml = optimize_eta(λ, µ, data);
    model = SSEconstant(λ, μ, ηml);

    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    rates = tree_rates(data, model, Fs, Ss);
    N = state_shifts(model, data, Ds, Ss; ape_order = false);
    nshift = sum(N, dims = (2,3))[:,1,1];
    append!(nshift, 0.0)
    rates[!,"nshift"] = nshift


    bf = posterior_prior_shift_odds(model,data)
    append!(bf, NaN)
    rates[!,"shift_bf"] = bf
    rates[!,"shift_bf_log"] = log10.(bf)

    #tip_rates(model, data, Ds, Fs)


    ## save data
    fpath = string("output/simulations/rate_shift_type/newick/", i, ".tre")
    writenewick(fpath, data, rates)

    fpath = string("output/simulations/rate_shift_type/rates/", i, ".csv")
    CSV.write(fpath, rates)

    fpath = string("output/simulations/rate_shift_type/jld2/", i, ".jld2")
    Nsum = sum(N, dims = 1)[1,:,:]
    save(fpath, 
        "N", N,
        "Nsum", Nsum,
        "lambda", λ,
        "mu", μ,
        "muml", μml,
        "lambdaml", λml,
        "etaml", ηml)
    next!(prog)
end
finish!(prog)



