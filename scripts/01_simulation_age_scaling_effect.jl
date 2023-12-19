using Revise
using BirthDeathSimulation
using ProgressMeter
using JLD2
using Distributions
using Random

## simulate trees conditional on survival until time t
Random.seed!(1234)

## a model with moderate rate variation
r0 = 0.04
r = [r0, 0.07, 0.10]
ϵ = 2/3
λ = r ./ (1 - ϵ)
μ = λ .- r
η = r0 / 50

model = bdsmodel(λ, μ, η)

maxtaxa = 50000
starting = [1]
n_trees = 500
n_states = length(model.λ)

tree_heights = collect(range(30, 100; length = 8))
n_heights = length(tree_heights)

trees = Array{Tree, 2}(undef, n_trees, n_heights)
N = Array{Int64, 4}(undef, n_trees, n_heights, n_states, n_states)

for j in 1:n_heights
    h = tree_heights[j]
    prog = ProgressUnknown("Complete trees sim (tree height = $h):")
    i = 1
    while i <= n_trees
        maxtime = tree_heights[j]
        state = 1
        tree = sim_bdshift(model, maxtime, maxtaxa, state)
        if length(tree.Leaves) < maxtaxa ## reject complete trees that termined when too many taxa
            if length(tree.Leaves) > 5 ## Reject completete trees where all taxa went extinct
                prune_extinct!(tree)
                if abs(treeheight(tree) - maxtime) < 0.001 ## reject complete trees where both root children did not survive
                    N0 = +([branch.N for (idx, branch) in tree.Branches]...)


                    if sum(N0) > 0 ## reject trees without any shifts
                        N[i,j,:,:] .= N0
                        trees[i,j] = tree
                        i += 1
                     end
                end
            end
        end
        ProgressMeter.next!(prog)
    end
    ProgressMeter.finish!(prog)
    sleep(0.1)
end


Ns = sum(N, dims = (3,4))[:,:,1,1]


prog = Progress(n_heights * n_trees, "Writing trees:")
for i in 1:n_heights
    h = string(Int64(tree_heights[i]))
    for j in 1:n_trees
        fpath = string(
            "data/simulations/age_scaling_effect/h", 
            h, "_", string(j), ".tre")
        writenewick(fpath, trees[j,i], model)
        next!(prog)
    end
end
finish!(prog)


using StatsPlots
x = [ntaxa(tree) for tree in trees]
histogram(vcat(x...))
y = hcat([tree_heights for i in 1:n_trees]...)
scatter(y', x, xlabel = "tree height", ylabel = "ntaxa")





