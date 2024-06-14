using Revise
using Pesto
using JLD2
using ProgressMeter
using Glob
using DataFrames
using CSV

fpaths = Glob.glob("output/empirical_joint/newick/*.tre")

inference = "empirical_joint"

df = CSV.read("output/empirical_munged.csv", DataFrame)
df = df[df[!,:inference] .== inference,:]

d = Dict()
names = unique(df[!,:name])
fpaths = Glob.glob("output/" * inference * "/jld2/*.jld2")
models = Dict{String, SSEconstant}()
@showprogress for fpath in fpaths
    name = split(Base.basename(fpath), ".")[1]

    λ = JLD2.load(fpath, "lambda")
    μ = JLD2.load(fpath, "mu")
    η = JLD2.load(fpath, "etaml")
    model = SSEconstant(λ, μ, η)
    models[name] = model
end


fpaths = ["data/empirical/" * name * ".tree" for name in names]
df = CSV.read("data/empirical/Phylogenies for DeepILS project.csv", DataFrame)
ρs = Dict()
for row in eachrow(df)
    fn = row["Filename"]

    if !ismissing(fn)
        name = replace(fn, ".tree" => "")
        ρ = row["P extant sampling"]
        ρs[name] = ρ
    end
end

trees = Dict()
datasets = Dict()
for fpath in fpaths
    println(fpath)
    tree = readtree(fpath)
    
    if all(tree.edge_length .> 0)
        bn = Base.Filesystem.basename(fpath)
        name = replace(bn, ".tree" => "")
        trees[name] = tree
        datasets[name] = SSEdata(trees[name], ρs[name])
    end
end

times = Dict{String,Vector{Float64}}()
mean_shift_rate = Dict{String,Vector{Float64}}()

@showprogress for name in keys(models)
    x, sr = Pesto.shift_rate_through_time(models[name], datasets[name])
    times[name] = x
    mean_shift_rate[name] = sr
end


for name in keys(models)
    df = DataFrame(
        "times" => times[name], 
        "nshift" => mean_shift_rate[name]) 
    fpath = "output/empirical_joint/shift_rate_through_time/$name.csv"
    CSV.write(fpath, df) 
end

