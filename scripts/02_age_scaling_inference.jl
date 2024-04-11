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
    * Branch rates CSV file
    * Estimate for η
=#

tree_heights = Int.(collect(range(30, 100; length = 8)))
n_iters = 500
n_heights = length(tree_heights)

io = open("output/prog_age_scaling_effect.txt", "w")
prog = ProgressMeter.Progress(n_iters * n_heights; desc = "Inference (age scaling effect): ", output = io);

for i in 1:n_iters
    for height in tree_heights[1:end]
        phy = readtree(string("data/simulations/age_scaling_effect/h", height, "_", i ,".tre"))
        ρ = 1.0
        data = SSEdata(phy, ρ);

        if length(data.tiplab) < 50 ## for small trees none of this works really
            continue
        end

        ## skip if already did this one
        fpath_jld2 = string("output/simulations/age_scaling_effect/jld2/h", height, "_", i, ".jld2")
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
            fpath = string("output/simulations/age_scaling_effect/newick/h", height, "_", i ,".tre")
            writenewick(fpath, data, rates)

            fpath = string("output/simulations/age_scaling_effect/rates/h", height, "_", i ,".csv")
            CSV.write(fpath, rates)

            fpath = string("output/simulations/age_scaling_effect/jld2/h", height, "_", i ,".jld2")
            Nsum = sum(N, dims = 1)[1,:,:]
            save(fpath, 
                "N", N,
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
end
finish!(prog)



