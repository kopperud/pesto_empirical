using Revise
using Pesto
using Glob
using ProgressMeter
using CSV
using Distributions
using JLD2
using DataFrames


fpaths = Glob.glob("data/simulations/up_vs_down/*.tre")

datasets = Dict{Int64, SSEdata}()

## read trees
@showprogress for fpath in fpaths
    ρ = 1.0
    tree = readtree(fpath)
    tree_index = parse(Int64, split(Base.basename(fpath), ".")[1])
    data = SSEdata(tree, ρ)
    datasets[tree_index] = data
end

## do the inference
n_iters = length(datasets)
io = open("output/prog_up_vs_down_joint2.txt","w")
prog = ProgressMeter.Progress(n_iters; desc = "Inference (joint r): ", output = io);
for (i, data) in datasets

    if length(data.tiplab) < 5 ## for small trees none of this works really
        continue
    end

    ## if already did this one
    fpath_jld2 = string("output/simulations/up_vs_down_joint2/jld2/", i, ".jld2")
    if isfile(fpath_jld2)
        continue
    end

    upper = [0.4, 2.0, 1.0]
    optres, model, n_attempts = optimize_hyperparameters(data; upper = upper)

    g,h = logistic(upper, 0.5)
    λml = g(sum(optres.minimizer))
    μml = g(sum(optres.minimizer[1:2]))

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
    fpath = string("output/simulations/up_vs_down_joint2/newick/", i, ".tre")
    writenewick(fpath, data, rates)

    fpath = string("output/simulations/up_vs_down_joint2/rates/", i, ".csv")
    CSV.write(fpath, rates)

    fpath = string("output/simulations/up_vs_down_joint2/jld2/", i, ".jld2")
    Nsum = sum(N, dims = 1)[1,:,:]
    save(fpath, 
        "N", N,
        "λml", λml,
        "μml", μml,
        "n_attempts", n_attempts,
        "λ", λ,
        "μ", μ,
        "ntip", ntip,
        "Nsum", Nsum,
        "etaml", ηml)
    next!(prog)
end
finish!(prog)



#=
fpaths2 = Glob.glob("output/simulations/up_vs_down_jointhyper/newick/*.tre")
for (i, data) in datasets

    fpath = string("output/simulations/up_vs_down_jointhyper/newick/", i, ".tre")

    ntip = length(data.tiplab)

    if ntip >= 5
        if !(fpath ∈ fpaths2)
            println(fpath)
            break
        end
    end

end
=#

close(io)

exit()
