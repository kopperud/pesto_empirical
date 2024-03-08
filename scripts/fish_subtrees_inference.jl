## Inference of the simulated datasets

using Revise
using Pesto
using Distributions
using CSV
using JLD2
using ProgressMeter
using Glob

#= 
These are the outputs I want for the inference:
    * Newick string
    * N matrix (3-dim)
    * Tip rates ??
    * Branch rates CSV file
    * Constant-rate Birth-Death Process estimates
    * Estimate for η
=#

fpaths = Glob.glob("data/empirical_fish_subtrees/*.tre")

prog = ProgressMeter.Progress(length(fpaths); desc = "Inference: ");
for i in eachindex(fpaths)
    phy = readtree(fpaths[i])
    node = parse(Int64,split(replace(Base.basename(fpaths[i]), ".tre" => ""), "_")[end])
    ntip = length(phy.tip_label)

    ρ = 1.0
    data = SSEdata(phy, ρ);

    λml, μml = estimate_constant_bdp(data)

    H = 0.587
    n = 7

    dλ = LogNormal(log(λml), H)
    dμ = LogNormal(log(µml), H)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)

    ηml = optimize_eta(λ, µ, data);
    treelength = sum(data.branch_lengths);
    model = SSEconstant(λ, μ, ηml);

    #Ds, Fs = backwards_forwards_pass(model, data);
    #Ss = ancestral_state_probabilities(data, Ds, Fs);

    #rates = tree_rates(data, model, Fs, Ss);
    rates = birth_death_shift(model, data)
    #N = state_shifts(model, data, Ds, Fs);

    ## save data
    fpath = string("output/empirical_fish_subtrees/newick/Actinopterygii_node_", node ,".tre")
    writenewick(fpath, data, rates)

    fpath = string("output/empirical_fish_subtrees/rates/Actinopterygii_node_", node ,".csv")
    CSV.write(fpath, rates)

    fpath = string("output/empirical_fish_subtrees/jld2/Actinopterygii_node_", node ,".jld2")
    #Nsum = sum(N, dims = 1)[1,:,:]
    treeheight = maximum(phy.node_depths)
    save(fpath, 
        #"N", N,
        #"Nsum", Nsum,
        "lambda", λ,
        "mu", μ,
        "muml", μml,
        "lambdaml", λml,
        "etaml", ηml,
        "ntip", ntip,
        "treelength", treelength,
        "treeheight", treeheight)
    next!(prog)
end
finish!(prog)



