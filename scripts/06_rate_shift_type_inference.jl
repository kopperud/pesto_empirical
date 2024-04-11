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
io = open("output/prog_rate_shift_type.txt", "w")
prog = ProgressMeter.Progress(n_iters; desc = "Inference (rate shift type): ", output = io);
for i in 1:n_iters
    phy = readtree(string("data/simulations/rate_shift_type/", i ,".tre"))
    ρ = 1.0
    data = SSEdata(phy, ρ);
    height = maximum(data.node_depth)

    if length(data.tiplab) < 50 ## for small trees none of this works really
        continue
    end

    ## skip if already did this one
    fpath_jld2 = string("output/simulations/rate_shift_type/jld2/", i, ".jld2")
    if isfile(fpath_jld2)
        continue
    end


    lower = [1e-08, 1e-04, 1e-04]
    upper = [0.3, 1.0, 1.0]


    try
        optres, model, n_attempts = optimize_hyperparameters(data; upper = upper, n_attempts = 10)

        g,h = logistic(lower, upper, 0.5)

        x = g(optres.minimizer)
        μml = sum(x[2])
        λml = maximum([5*x[1], x[2]]) + x[3]

        ntip = length(data.tiplab)
        treelength = sum(data.branch_lengths);

        λ = model.λ
        μ = model.μ
        ηml = model.η

        Ds, Fs = backwards_forwards_pass(model, data);
        Ss = ancestral_state_probabilities(data, Ds, Fs);

        rates = tree_rates(data, model, Fs, Ss);
        N = state_shifts(model, data, Ds, Ss);
        nshift = sum(N, dims = (2,3))[:,1,1];
        append!(nshift, 0.0)
        rates[!,"nshift"] = nshift


        bf = posterior_prior_shift_odds(model,data)
        append!(bf, NaN)
        rates[!,"shift_bf"] = bf
        rates[!,"shift_bf_log"] = log10.(bf)

 
        ## save data
        fpath = string("output/simulations/rate_shift_type/newick/", i, ".tre")
        writenewick(fpath, data, rates)

        fpath = string("output/simulations/rate_shift_type/rates/", i, ".csv")
        CSV.write(fpath, rates)

        fpath = string("output/simulations/rate_shift_type/jld2/", i, ".jld2")
        Nsum = sum(N, dims = 1)[1,:,:]
        save(fpath, 
            ##"N", N, ## do I need the full N matrix?
            "Nsum", Nsum,
            "lambda", λ,
            "mu", μ,
            "muml", μml,
            "lambdaml", λml,
            "n_attempts", n_attempts,
            "ntip", ntip,
            "etaml", ηml,
            "treeheight", height,
            "treelength", treelength)
      
    catch e
        if e isa Pesto.ConvergenceException
            continue
        else
            rethrow(e)
        end
    end

   next!(prog)
end
finish!(prog)



