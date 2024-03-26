using Revise
using Pesto
using Glob
using ProgressMeter
using CSV
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
io = open("output/prog_up_vs_down.txt","w")
prog = ProgressMeter.Progress(n_iters; desc = "Inference (known r): ", output = io);
for (i, data) in datasets

    if length(data.tiplab) < 5 ## for small trees none of this works really
        continue
    end

    ## if already did this one
    fpath_jld2 = string("output/simulations/up_vs_down/jld2/", i, ".jld2")
    if isfile(fpath_jld2)
        continue
    end

    r = [0.01, 0.03, 0.05, 0.07, 0.09]
    ϵ = 2/3
    λ = r ./ (1-ϵ)
    μ = λ .- r

    ηml = optimize_eta(λ, µ, data);
    model = SSEconstant(λ, μ, ηml);

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
    fpath = string("output/simulations/up_vs_down/newick/", i, ".tre")
    writenewick(fpath, data, rates)

    fpath = string("output/simulations/up_vs_down/rates/", i, ".csv")
    CSV.write(fpath, rates)

    fpath = string("output/simulations/up_vs_down/jld2/", i, ".jld2")
    Nsum = sum(N, dims = 1)[1,:,:]
    save(fpath, 
        "N", N,
        "ntip", ntip,
        "Nsum", Nsum,
        "etaml", ηml)
    next!(prog)
end
finish!(prog)


close(io)

exit()

