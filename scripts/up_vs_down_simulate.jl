using Revise
using BirthDeathSimulation
using LinearAlgebra
using CairoMakie
using Printf
using LaTeXStrings
import Random

Random.seed!(123)

## five state model
r = [0.01, 0.03, 0.05, 0.07, 0.09]
#r = [0.02, 0.02, 0.02, 0.02, 0.02]
Δr = r .- r'
ϵ = 2/3
λ = r ./ (1-ϵ)
μ = λ .- r

η = 0.001
model = bdsmodel(λ, µ, η)

max_time = 100.0
max_taxa = 10_000
starting_state = 3

trees = Tree[]

ntrees = 0
while ntrees < 1500
    #starting_state = sample(1:5,1)[1]
    tree = sim_bdshift(model, max_time, max_taxa, starting_state)
    #if ntaxa(tree) > 1 ## condition on that the tree left atleast 2 living taxa
    if ntaxa(tree) < max_taxa
        ntrees += 1
        push!(trees, tree)
    end
end


nt = [ntaxa(tree) for tree in trees]

fig = Figure()
ax1 = Axis(fig[1,1], 
    ylabel = "number of trees", 
    title = "a) extinct trees", 
    xticks = (Float64[0], String["0"]),
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    )
ax2 = Axis(fig[1,2], #xlabel = "number of taxa", 
    ylabel = "", title = "b) trees that did not go exintct",
        xticks = (0:4, ["1", "10", "100", "1000", "10000"]),
        yticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
        )
xlabel = Label(fig[2,1:2], "number of taxa")
CairoMakie.barplot!(ax1, [0], [sum(nt .== 0)], width = 0.5, color = :gray)
CairoMakie.hist!(ax2, log10.([x for x in nt if x > 0]), color = :gray)
linkyaxes!(ax1, ax2)
colsize!(fig.layout, 1, Relative(0.1))
colsize!(fig.layout, 2, Relative(0.9))
fig


function rate_shift_matrix(tree::Tree, model)
    K = length(model.λ)
    N = zeros(Int64, K, K)

    for (key, branch) in tree.Branches
        N[:,:] += branch.N
    end 
    return(N)
end

K = length(model.λ)

Ncomplete = zeros(Int64, ntrees, K, K)
Nreconstructed = zeros(Int64, ntrees, K, K)
magnitudes_complete = Float64[]
magnitudes_reconstructed = Float64[]

reconstructed_trees = Tree[]

for (i, tree) in enumerate(trees)
    N = rate_shift_matrix(trees[i], model)
    Ncomplete[i,:,:] .= N

    if sum(N) > 0
        mag = sum(N .* Δr) / sum(N)
        append!(magnitudes_complete, mag)
    end

    if ntaxa(trees[i]) > 1 ## can't do it for just 1 surviving lineage, not supported
        ## can argue that 1 lineage is not a phylogeny anyway
        ## define phylogeny := at least two lineages that had a common ancestor
        prune_extinct!(trees[i])
        push!(reconstructed_trees, trees[i])
        N = rate_shift_matrix(trees[i], model)
        Nreconstructed[i,:,:] .= N

        if sum(N) > 0
            mag = sum(N .* Δr) / sum(N)
            append!(magnitudes_reconstructed, mag)
        end
    end
end




## save trees to file
for (i, tree) in enumerate(reconstructed_trees)
    fpath = string("data/simulations/up_vs_down/",i,".tre")
    writenewick(fpath, tree, model)
end

## tree height are not equal, because either the left or right subtree went extinct,
## and then the branch leading up to the MRCA was pruned
treeheights = [treeheight(tree) for tree in reconstructed_trees]
hist(treeheights)




Nc = sum(Ncomplete, dims = 1)[1,:,:]
Nr = sum(Nreconstructed, dims = 1)[1,:,:]


for i in 1:ntrees
    sum(N .* Δr) ./ sum(N)
end


println("Number of shifts for complete phylogeny:")
println("increasing net-div: \t", sum(tril(Nc)))
println("decreasing net-div: \t", sum(triu(Nc)))
println("\n")
println("Number of shifts for reconstructed phylogeny:")
println("increasing net-div: \t", sum(tril(Nr)))
println("decreasing net-div: \t", sum(triu(Nr)))
println("\n")



fig2 = Figure(size = (1000, 600))

Δr_v = round.(vcat(Δr...), digits = 3)

xt = [-0.08, -0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06, 0.08]
xtl = [@sprintf "%2.2f" x for x in xt]

ax1 = Axis(
    fig2[1,1], xticks = (xt, xtl),
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    ylabel = L"\text{number of rate shifts }(N)",
    xlabel = L"\text{shift size in net diversification rate }(\Delta r)",
    )
ax2 = Axis(
    fig2[1,2], xticks = (xt, xtl),
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    ylabel = L"\text{number of rate shifts }(N)",
    xlabel = L"\text{shift size in net diversification rate }(\Delta r)",
    )

y = vcat(Nc...)
barplot!(ax1, Δr_v, y, label = "")

y = vcat(Nr...)
barplot!(ax2, Δr_v, y, label = "")

linkaxes!(ax1, ax2)
fig2

#########################
##
##    magnitudes
##
#########################

fig3 = Figure(size = (500, 600))

ax1 = Axis(
    fig3[1,1], #xticks = (xt, xtl),
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    ylabel = L"\text{number of trees}",
    title = "complete trees (≥ 1 shift)",
    )
ax2 = Axis(
    fig3[2,1], #xticks = (xt, xtl),
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    title = "reconstructed trees (≥ 1 shift)",
    ylabel = L"\text{number of trees}",
    xlabel = L"\text{magnitude }=\frac{\sum_{i,j}(r_i-r_j)N_{ij}}{\sum_{i,j}N_{ij}}",
    )

hist!(ax1, magnitudes_complete, color = :gray)
lines!(ax1, [0.0, 0.0], [0.0, 150], color = :red, linestyle = :dash)
hist!(ax2, magnitudes_reconstructed, color = :gray)
lines!(ax2, [0.0, 0.0], [0.0, 150], color = :red, linestyle = :dash)
linkaxes!(ax1, ax2)
fig3

















#= 
P(+|i=1) = 1
P(+|i=2) = 3/4
P(+|i=3) = 1/2
P(+|i=4) = 1/4
P(+|i=5) = 0 
 =#















Ncomplete
Nsurvived

Δr = r .- r'

## barplot





