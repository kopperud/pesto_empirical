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

tree_heights = collect(range(30, 100; length = 8))
n_iters = 500
n_heights = length(tree_heights)

prog = ProgressMeter.Progress(n_iters * n_heights; desc = "Inference: ");
for i in 1:n_iters
    for h in tree_heights[1:end]
        phy = readtree(string("data/simulations/age_scaling_effect/h", Int64(h), "_", i ,".tre"))
        ρ = 1.0
        data = SSEdata(phy, ρ);

        upper = [0.4, 2.0, 1.0]

        try
            optres, model, n_attempts = optimize_hyperparameters(data; upper = upper, n_attempts = 10)

            g,h = logistic(upper, 0.5)

            x = g(optres.minimizer)
            μml = sum(x[1:2])
            λml = sum(x)

            ntip = length(data.tiplab)

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
            fpath = string("output/simulations/age_scaling_effect/newick/h", Int64(h), "_", i ,".tre")
            writenewick(fpath, data, rates)

            fpath = string("output/simulations/age_scaling_effect/rates/h", Int64(h), "_", i ,".csv")
            CSV.write(fpath, rates)

            fpath = string("output/simulations/age_scaling_effect/jld2/h", Int64(h), "_", i ,".jld2")
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



