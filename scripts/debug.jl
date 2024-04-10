using Revise
using Pesto
using ProgressMeter
using Glob

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


for (i, data) in datasets

    if length(data.tiplab) < 5 ## for small trees none of this works really
        continue
    end

    ## if already did this one
    fpath_jld2 = string("output/simulations/up_vs_down_joint2/jld2/", i, ".jld2")
    if isfile(fpath_jld2)
        continue
    else
        println(fpath_jld2)
        break
    end
end


upper = [0.4, 2.0, 1.0]
optres, model, n_attempts = optimize_hyperparameters(datasets[716]; upper = upper, n_attempts = 3)

optres, model, n_attempts = optimize_hyperparameters(datasets[823]; upper = upper, n_attempts = 3)


d = Dict(i => (length(x.tiplab)-20)^2 for (i,x) in datasets)

Dict(i => n for (i, n) in d if n < 50)

d = Dict(i => length(x.tiplab) for (i, x) in datasets)
m = Dict(i => n for (i,n) in d if (n > 20) & (n < 50))



optres, model, n_attempts = optimize_hyperparameters(datasets[278]; upper = upper, n_attempts = 3)









