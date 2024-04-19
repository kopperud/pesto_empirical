using Revise
using Glob
using ProgressMeter
using CSV
using JLD2
using DataFrames
using RCall
using CairoMakie
using Statistics


function rm_na(x)
    res = x[.!isnan.(x)]
    return(res)
end


function find_tree_index(fpath::String)
    tree_index = parse(Int64, split(basename(fpath), ".")[1])
    return(tree_index)
end

datapaths_joint2 = Glob.glob("output/simulations/up_vs_down_joint2/jld2/*.jld2")
treepaths = Glob.glob("data/simulations/up_vs_down/*.tre")
ntrees = length(treepaths)


tree_lengths = Dict{Int64, Float64}()
tree_heights = zeros(ntrees)
tree_heights[:] .= NaN

@showprogress for i in eachindex(treepaths)
    fpath = treepaths[i]
    tree_index = parse(Int64, split(basename(fpath), ".")[1])
    @rput fpath
    R"""
    library(ape)
    tr <- read.tree(fpath)
    tl <- sum(tr$edge.length)
    th <- max(node.depth.edgelength(tr))
    """
    @rget tl
    @rget th

    tree_heights[tree_index] = th
    tree_lengths[tree_index] = tl
end

@rput treepaths
R""" ## this is kind of slow and not necessary, I didn't save the N matrices before
library(treeio)

items <- list()
ntrees <- length(treepaths)

for (i in 1:ntrees){
    tree_index <- as.numeric(strsplit(basename(treepaths[i]), "\\.")[[1]][1])

    tr <- read.beast.newick(treepaths[i])
    l <- lapply(tr@data$N, function(x) matrix(x, 5, 5))
    N <- Reduce("+", l)
    items[[tree_index]] <- N
}
"""
@rget items

r = [0.01, 0.03, 0.05, 0.07, 0.09]
Δr = r .- r'
Ntrue = zeros(Int64, ntrees, 5, 5)
magnitudes_true = Float64[]
for tree_index in 1:ntrees
    Ntrue[tree_index,:,:] .= items[tree_index]
    mag = sum(Ntrue[tree_index,:,:] .* Δr) / sum(Ntrue[tree_index,:,:])
    push!(magnitudes_true, mag)
end


N_estimated_joint2 = zeros(Float64, ntrees, 36, 36)
N_estimated_joint2[:,:,:] .= NaN
magnitudes_estimated_joint2 = zeros(ntrees)
magnitudes_estimated_joint2[:] .= NaN
ntip = zeros(ntrees)
ntip[:] .= NaN
for (i, path) in enumerate(datapaths_joint2)
    tree_index = parse(Int64, split(basename(path), ".")[1])
    N_estimated_joint2[tree_index,:,:] .= load(path, "Nsum")

    λ = load(path, "λ")
    μ = load(path, "μ")
    r = λ .- μ
    Δr = r .- r'

    ntip[tree_index] = load(path, "ntip")

    mag = sum(N_estimated_joint2[tree_index,:,:] .* Δr) / sum(N_estimated_joint2[tree_index,:,:])
    magnitudes_estimated_joint2[tree_index] = mag
end


println("mean magnitude (true): \t $(mean(rm_na(magnitudes_true)))");
println("mean magnitude (estimated, unknown λ,μ,η, joint2): \t $(mean(rm_na(magnitudes_estimated_joint2)))")

println("median magnitude (true): \t $(median(rm_na(magnitudes_true)))");
println("median magnitude (estimated, unknown λ,μ,η, joint2): \t $(median(rm_na(magnitudes_estimated_joint2)))")



#################################
#
#      plot estimation error
#
############################


fig = Figure(size = (500, 400))

ax1 = Axis(fig[1,1], 
        xlabel = L"\text{magnitude}",
        title = L"\text{a) true}",
        xgridvisible = false, ygridvisible = false,
        topspinevisible = false, rightspinevisible = false,
        ylabel = L"\text{number of trees}")

ax2 = Axis(fig[1,2], 
        xlabel = L"\text{magnitude}",
        title = L"\text{b) estimated}",
        topspinevisible = false, rightspinevisible = false,
        xgridvisible = false, ygridvisible = false)

hist!(ax1, rm_na(magnitudes_true), bins = 15, strokecolor = :black, strokewidth = 1)
hist!(ax2, rm_na(magnitudes_estimated_joint2), bins = 30, strokecolor = :black, strokewidth = 1)
linkaxes!(ax1, ax2)

error = magnitudes_estimated_joint2 .- magnitudes_true

df = DataFrame(
    "ntip" => ntip,
    "error" => error,
    "height" => tree_heights,
    "magnitude_estimated" => magnitudes_estimated_joint2,
    "magnitude_true" => magnitudes_true,
)
df = filter(:error => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)


ax3 = Axis(fig[2,1], xgridvisible = false, ygridvisible = false, 
        topspinevisible = false, rightspinevisible = false,
        ylabel = L"\text{number of trees}",
        title = L"\text{c) error distribution}",
        xlabel = L"\text{estimation error (magnitude)}")

hist!(ax3, df[!,:error], bins = 30, color = :orange, strokecolor = :black, strokewidth = 1, label = "error")
lines!(ax3, repeat([0.0], 2), [0.0, 150], linestyle = :dash, color = :red, label = "zero error")

xt = [50, 500, 5000]
ax4 = Axis(fig[2,2],
    xgridvisible = false, ygridvisible = false, 
    topspinevisible = false, rightspinevisible = false,
    xticks = xt,
    xscale = log10,
    title = L"\text{d) error vs tips}",
    xlabel = L"\text{number of tips}",
    ylabel = L"\text{estimation error (magnitude)}")
    
scatter!(ax4, df[!,:ntip], df[!,:error], markersize = 5, color = (:orange, 0.5))
lines!(ax4, [extrema(df[!,:ntip])...], repeat([0.0], 2), linestyle = :dash, color = :red, label = "zero error")

fig
save("figures/magnitude_results.pdf", fig)





fig2 = Figure(size = (450, 350))
ax1 = Axis(fig2[1,1],
    rightspinevisible = false,
    topspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    #xlabel = L"\text{tree height (Ma)}",
    ylabel = L"\text{magnitude}",
    title = L"\text{a) true}")

    scatter!(ax1, df[!,:height], df[!,:magnitude_true], markersize = 5)

ax2 = Axis(fig2[1,2],
    rightspinevisible = false,
    topspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    #xlabel = L"\text{tree height (Ma)}",
    #ylabel = L"\text{magnitude}",
    title = L"\text{b) estimated}")

linkaxes!(ax1, ax2)
scatter!(ax2, df[!,:height], df[!,:magnitude_estimated], markersize = 5)

ax3 = Axis(fig2[2,1:2],
    rightspinevisible = false,
    topspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    xlabel = L"\text{tree height (Ma)}",
    ylabel = L"\text{error (magnitude)}",
    title = L"\text{c) estimation error } = \text{mag}_\text{estimated} - \text{mag}_\text{true}")

scatter!(ax3, df[!,:height], df[!,:error], markersize = 5, color = (:orange, 0.5))
lines!(ax3, [extrema(df[!,:height])...], [0.0, 0.0], markersize = 5, color = :red, linestyle = :dash)

fig2
save("figures/magnitude_vs_height.pdf", fig2)


hist(df[!,:height])



