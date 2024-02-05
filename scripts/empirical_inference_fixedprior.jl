## Inference of the simulated datasets

using Revise
using Pesto
using Distributions
using CSV
using JLD2
using ProgressMeter
using Glob
using DataFrames

#= 
These are the outputs I want for the inference:
    * Newick string
    * N matrix (3-dim)
    * Tip rates ??
    * Branch rates CSV file
    * Constant-rate Birth-Death Process estimates
    * Estimate for η
=#



df = CSV.read("data/empirical/Phylogenies for DeepILS project.csv", DataFrame)
ρs = Dict{String, Float64}()
for row in eachrow(df)
    fn = row["Filename"]
    ρ = row["P extant sampling"]
    if ρ isa Real && fn isa String
        ρs[fn] = ρ
    end
end

fpaths = Glob.glob("data/empirical/*.tre*")
fpaths = [
    fpath for fpath in fpaths if Base.basename(fpath) ∈ keys(ρs)
]

skipped_trees = [
    "Seed_plants_Smith&Brown2018_GBMB", ## skip because not binary
    "Seed_plants_Smith&Brown2018_GBMB_binary", ## same, it's binary but shitty
    "Mimosoideae_Ringelberg2023", ## skip because negative branch lengths
    "Nymphalidae_Chazot2021", ## negative branch lenghts
    "Solanaceae_Sarkinen2013" ## ??
]

## verify that all trees can be read
for fpath in fpaths
    name = join(split(Base.basename(fpath), ".")[1:end-1])
    if name in skipped_trees
        continue
    end

    println(name)
    phy = readtree(fpath)
    data = SSEdata(phy, 1.0)
    println(fpath, ",\t n = ", length(data.tiplab))
end

completed_jobs = [
    split(Base.basename(x), ".")[1] for x in Glob.glob("output/empirical_fixedprior/jld2/*.jld2")
]


n_iters = length(fpaths)
n_iters - length(completed_jobs)
prog = ProgressMeter.Progress(n_iters; desc = "Inference: ");
for fpath in fpaths
    phy = readtree(fpath)
    name = join(split(Base.basename(fpath), ".")[1:end-1])

    if name in skipped_trees
        continue
    end

    if name in completed_jobs
        continue
    end
    println(name)

    ## read in ρ
    ρ = ρs[basename(fpath)]
    data = SSEdata(phy, ρ);

    #λml, μml = estimate_constant_bdp(data)
    λml = 0.15
    μml = 0.10

    H = 0.587
    n = 10

    dλ = LogNormal(log(λml), 2*H)
    dμ = LogNormal(log(µml), 2*H)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)

    ηml = optimize_eta(λ, µ, data);
    model = SSEconstant(λ, μ, ηml);
    logl = logL_root(model, data)

    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    rates = tree_rates(data, model, Fs, Ss);
    N = state_shifts(model, data, Ds, Fs);
    nshift = sum(N, dims = (2,3))[:,1,1];
    append!(nshift, 0.0)
    rates[!,"nshift"] = nshift


    bf = posterior_prior_shift_odds(model,data)
    append!(bf, NaN)
    rates[!,"shift_bf"] = bf
    rates[!,"shift_bf_log"] = log10.(bf)

    #tip_rates(model, data, Ds, Fs)


    ## save data
    fpath = string("output/empirical_fixedprior/newick/", name, ".tre")
    writenewick(fpath, data, rates)

    fpath = string("output/empirical_fixedprior/rates/", name, ".csv")
    CSV.write(fpath, rates)

    fpath = string("output/empirical_fixedprior/jld2/", name, ".jld2")
    Nsum = sum(N, dims = 1)[1,:,:]
    save(fpath, 
        "N", N,
        "Nsum", Nsum,
        "lambda", λ,
        "logl", logl,
        "mu", μ,
        "muml", μml,
        "lambdaml", λml,
        "etaml", ηml)
    next!(prog)
end
finish!(prog)




