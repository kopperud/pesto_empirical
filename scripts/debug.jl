using Revise
using Pesto
using ProgressMeter
using Glob

fpaths = Glob.glob("data/simulations/age_scaling_effect/*.tre")

datasets = Dict{Tuple{Int64, Int64}, SSEdata}()

## read trees
@showprogress for fpath in fpaths
    ρ = 1.0
    tree = readtree(fpath)
    bn = split(Base.basename(fpath), ".")[1]
    a, b = split(bn, "_")
    height = parse(Int64, replace(a, "h" => ""))
    tree_index = parse(Int64, b)
    #height = 
    #tree_index = parse(Int64, split(Base.basename(fpath), ".")[1])
    data = SSEdata(tree, ρ)
    datasets[height, tree_index] = data
end

model_eb = empirical_bayes(datasets[100,11])
logL_root(model_eb, datasets[100, 11])

models = SSEconstant[]
optres = []
is = []
logls = []

for i in 1:5
    or, model, i = optimize_hyperparameters(datasets[100,11])
    push!(optres, or)
    push!(models, model)
    push!(is, i)
    push!(logls, logl_root(model, datasets[100,11]))
end
    #logL_root(model, datasets[100,11])


rates = birth_death_shift(model, datasets[100,11])

using CairoMakie
treeplot(datasets[100,11], rates)


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









