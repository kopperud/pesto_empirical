using Revise
using Glob
using ProgressMeter
using CSV
using JLD2
using DataFrames
using RCall
using CairoMakie



datapaths_known = Glob.glob("output/simulations/up_vs_down/jld2/*.jld2")

datapaths_twostep = Glob.glob("output/simulations/up_vs_down_twostep/jld2/*.jld2")

datapaths_joint = Glob.glob("output/simulations/up_vs_down_joint/jld2/*.jld2")

treepaths = Glob.glob("data/simulations/up_vs_down/*.tre")


tree_lengths = Dict{Int64, Float64}()
@showprogress for i in eachindex(treepaths)
    fpath = treepaths[i]
    tree_index = parse(Int64, split(basename(fpath), ".")[1])
    @rput fpath
    R"""
    tr <- read.tree(fpath)
    tl <- sum(tr$edge.length)
    """
    @rget tl
    tree_lengths[tree_index] = tl
end

## true

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
ntrees = length(treepaths)
Ntrue = zeros(Int64, ntrees, 5, 5)
magnitudes_true = Float64[]
for tree_index in 1:ntrees
    Ntrue[tree_index,:,:] .= items[tree_index]
    mag = sum(Ntrue[tree_index,:,:] .* Δr) / sum(Ntrue[tree_index,:,:])
    #if !isnan(mag)
    push!(magnitudes_true, mag)
    #end
end

#Δr
r = [0.01, 0.03, 0.05, 0.07, 0.09]
Δr = r .- r'
N_estimated_known = zeros(Float64, ntrees, 5, 5)
magnitudes_estimated_known = zeros(ntrees)
magnitudes_estimated_known[:] .= NaN
for (i, path) in enumerate(datapaths_known)
    tree_index = parse(Int64, split(basename(path), ".")[1])
    N_estimated_known[tree_index,:,:] .= load(path, "Nsum")

    mag = sum(N_estimated_known[tree_index,:,:] .* Δr) / sum(N_estimated_known[tree_index,:,:])

    #push!(magnitudes_estimated_known, mag)
    magnitudes_estimated_known[tree_index] = mag
end

fig = Figure(size = (400, 600))
ax1 = Axis(fig[1,1], xscale = log10, xlabel = "number of tips", ylabel = "magnitude")
ax2 = Axis(fig[2,1], xscale = log10, xlabel = "number of tips", ylabel = "magnitude")
scatter!(ax1, ntip_all, magnitudes_estimated_known)
scatter!(ax2, ntip_all, magnitudes_true)
linkaxes!(ax1, ax2)

fig


Δrs_twostep = zeros(Float64, ntrees, 36, 36)
N_estimated_twostep = zeros(Float64, ntrees, 36, 36)
magnitudes_estimated_twostep = zeros(ntrees)
magnitudes_estimated_twostep[:] .= NaN
N_estimated_twostep[:,:,:] .= NaN
for (i, path) in enumerate(datapaths_twostep)
    tree_index = parse(Int64, split(basename(path), ".")[1])
    N_estimated_twostep[tree_index,:,:] .= load(path, "Nsum")

    λ = load(path, "λ")
    μ = load(path, "μ")
    r = λ .- μ
    Δr = r .- r'

    mag = sum(N_estimated_twostep[tree_index,:,:] .* Δr) / sum(N_estimated_twostep[tree_index,:,:])

    magnitudes_estimated_twostep[tree_index] = mag
end

N_estimated_joint = zeros(Float64, ntrees, 36, 36)
N_estimated_joint[:,:,:] .= NaN
magnitudes_estimated_joint = zeros(ntrees)
magnitudes_estimated_joint[:] .= NaN
for (i, path) in enumerate(datapaths_joint)
    tree_index = parse(Int64, split(basename(path), ".")[1])
    N_estimated_joint[tree_index,:,:] .= load(path, "Nsum")

    λ = load(path, "λ")
    μ = load(path, "μ")
    r = λ .- μ
    Δr = r .- r'

    mag = sum(N_estimated_joint[tree_index,:,:] .* Δr) / sum(N_estimated_joint[tree_index,:,:])
    
    magnitudes_estimated_joint[tree_index] = mag
end



fig = Figure(size = (650, 450))

ax1 = Axis(fig[1,1], 
    title = L"\text{true shifts (reconstructed trees)}",
    xgridvisible = false, ygridvisible = false,
    ylabel = L"\text{number of trees}")
ax2 = Axis(fig[2,1], 
    title = L"\text{estimated shifts, known }𝛌,𝛍, \text{ unknown }η",
    xlabel = "magnitude",
    xgridvisible = false, ygridvisible = false,
    ylabel = L"\text{number of trees}")
ax3 = Axis(fig[1,2], 
    title = L"\text{estimated shifts, unknown }𝛌,𝛍,η\text{ (two step})",
    xgridvisible = false, ygridvisible = false,
    ylabel = L"\text{number of trees}")
ax4 = Axis(fig[2,2], 
        title = L"\text{estimated shifts, unknown }𝛌,𝛍,η\text{ (joint})",
        xlabel = "magnitude",
        xgridvisible = false, ygridvisible = false,
        ylabel = L"\text{number of trees}")
linkaxes!(ax1, ax2, ax3, ax4)
#linkxaxes!(ax1, ax4)
function rm_na(x)
    res = x[.!isnan.(x)]
    return(res)
end

hist!(ax1, rm_na(magnitudes_true), bins = 50, color = :black)
hist!(ax2, rm_na(magnitudes_estimated_known), bins = 50, color = :green)
hist!(ax3, rm_na(magnitudes_estimated_twostep), bins = 50, color = :orange)
hist!(ax4, rm_na(magnitudes_estimated_joint), bins = 50, color = :purple)
fig
#save("figures/up_vs_down_simulation.pdf", fig)
save("figures/Q2_simulations.pdf", fig)


function mean(x::Vector{T}) where {T <: Real}
    s = sum(x)
    n = length(x)
    m = s/n
    return(m)
end

function find_tree_index(fpath::String)
    tree_index = parse(Int64, split(basename(fpath), ".")[1])
    return(tree_index)
end

function find_number_of_tips(tree_index)
    fpath = string("data/simulations/up_vs_down/", tree_index, ".tre")
    @rput fpath
    R"""
    tr <- read.tree(fpath)
    ntip <- length(tr$tip.label)
    """
    @rget ntip
    return(ntip)
end





## eta estimates
η1 = [load(path, "etaml") for path in datapaths_known]
η2 = [load(path, "etaml") for path in datapaths_twostep]
η3 = [load(path, "etaml") for path in datapaths_joint]



#[load(path,"ntip") for path in datapaths_joint]
ntip_all = [find_number_of_tips(find_tree_index(fpath)) for fpath in Glob.glob("data/simulations/up_vs_down/*.tre")]
ntip1 = [find_number_of_tips(find_tree_index(fpath)) for fpath in datapaths_known]
ntip2 = [find_number_of_tips(find_tree_index(fpath)) for fpath in datapaths_twostep]
ntip3 = [find_number_of_tips(find_tree_index(fpath)) for fpath in datapaths_joint]

tree_lengths_all = []
for fpath in Glob.glob("data/simulations/up_vs_down/*.tre")
    #fpath = string("data/simulations/up_vs_down/", tree_index, ".tre")
    @rput fpath
    R"""
    tr <- read.tree(fpath)
    tl <- sum(tr$edge.length)
    """
    @rget tl
    push!(tree_lengths_all, tl)
end

println("mean magnitude (true): \t $(mean(magnitudes_true))");
println("mean magnitude (estimated, known λ,μ, unknown η): \t $(mean(magnitudes_estimated_known))")
println("mean magnitude (estimated, unknown λ,μ,η, two-step): \t $(mean(magnitudes_estimated_twostep))")
println("mean magnitude (estimated, unknown λ,μ,η, joint): \t $(mean(magnitudes_estimated_joint))")

using Statistics
println("median magnitude (true): \t $(median(magnitudes_true))");
println("median magnitude (estimated, known λ,μ, unknown η): \t $(median(magnitudes_estimated_known))")
println("median magnitude (estimated, unknown λ,μ,η, two-step): \t $(median(magnitudes_estimated_twostep))")
println("median magnitude (estimated, unknown λ,μ,η, joint): \t $(median(magnitudes_estimated_joint))")




fig2 = Figure(size = (600, 400));

ax1 = Axis(fig2[1,1], xscale = log10, yscale = log10, 
    title = L"\text{true}", xticklabelsvisible = false, 
    ylabel = "shift rate (η)",
    topspinevisible = false, rightspinevisible = false, 
    xgridvisible = false, ygridvisible = false)
ax2 = Axis(fig2[2,1], xscale = log10, yscale = log10, 
    title = L"\text{estimated }(𝛌,𝛍\text{ known}, η \text{ unknown)}", 
    ylabel = "shift rate (η)", topspinevisible = false, rightspinevisible = false, xgridvisible = false, ygridvisible = false,
    xlabel = "number of tips")
ax3 = Axis(fig2[1,2], xscale = log10, yscale = log10,
 title = L"\text{estimated }(𝛌,𝛍,η \text{ unknown}), \text{ two step}", 
 xticklabelsvisible = false, ylabel = "shift rate (η)", topspinevisible = false, rightspinevisible = false, xgridvisible = false, ygridvisible = false)
ax4 = Axis(fig2[2,2], xscale = log10, yscale = log10,
title = L"\text{estimated }(𝛌,𝛍,η \text{ unknown}), \text{ joint}",
 ylabel = "shift rate (η)", topspinevisible = false, rightspinevisible = false, xgridvisible = false, ygridvisible = false,
 xlabel = "number of tips")

linkaxes!(ax1, ax2, ax3, ax4)

true_shift_rate = [sum(item) for item in items] ./ tree_lengths_all
true_shift_rate[true_shift_rate .== 0] .= 1e-8

scatter!(ax1, ntip_all, true_shift_rate, color = :gray)
scatter!(ax2, ntip1, η1, bins = 50, color = :green)
scatter!(ax3, ntip2, η2, bins = 50, color = :orange)
scatter!(ax4, ntip3, η3, bins = 50, color = :purple)
for ax in (ax1, ax2, ax3, ax4)
    lines!(ax, [2.0, 10_000.0], [0.001, 0.001], color = :red, linestyle = :dash) 
end

fig2
save("figures/up_vs_down_shift_rate_estimates.pdf", fig2)


### number of attempts
attempts = [load(path, "n_attempts") for path in datapaths_joint]
λml = [load(path, "λml") for path in datapaths_joint]
μml = [load(path, "μml") for path in datapaths_joint]
rml = λml .- μml
ηml = [load(path, "etaml") for path in datapaths_joint]


fig3 = Figure(size = (700, 500))
ax1 = Axis(fig3[1,1], xgridvisible = false, ygridvisible = false,
    ylabel = "number of trees",
    xlabel = "attempts (convergence)")
hist!(ax1, attempts, bins = 10)

fig3

ax2 = Axis(fig3[1,2],
        xgridvisible = false, ygridvisible = false,
        xlabel = "number of tips",
        xscale = log10,
        ylabel = "attempts (convergence)")
scatter!(ax2, ntip3, attempts)

ax3 = Axis(fig3[2,1],
        xgridvisible = false, ygridvisible = false,
        xlabel = "r = λ - μ",
        xscale = log10,
        ylabel = "attempts (convergence)")
scatter!(ax3, rml, attempts)

ax4 = Axis(fig3[2,2],
        xgridvisible = false, ygridvisible = false,
        xlabel = "shift rate (η)",
        xscale = log10,
        ylabel = "attempts (convergence)")
scatter!(ax4, ηml, attempts)

ax5 = Axis(fig3[1,3],
        xgridvisible = false, ygridvisible = false,
        xlabel = "speciation (λ)",
        xscale = log10,
        ylabel = "attempts (convergence)")
scatter!(ax5, λml, attempts)

ax6 = Axis(fig3[2,3],
        xgridvisible = false, ygridvisible = false,
        xlabel = "extinction (μ)",
        xscale = log10,
        ylabel = "attempts (convergence)")
scatter!(ax6, μml, attempts)

save("figures/number_of_attempts_diagnostic.pdf", fig3)
fig3

## error


fig = Figure(size = (400, 550))
ax1 = Axis(fig[1,1], xgridvisible = false, ygridvisible = false, title = L"\text{known }𝛌,𝛍,\text{ unknown } η")
ax2 = Axis(fig[2,1], xgridvisible = false, ygridvisible = false, title = L"\text{unknown }𝛌,𝛍,η\text{, two step}")
ax3 = Axis(fig[3,1], xgridvisible = false, ygridvisible = false, title = L"\text{unknown }𝛌,𝛍,η\text{, joint}", 
        xlabel = L"\text{error} = \text{mag}_\text{true} - \text{mag}_\text{estimated}")

mags = (magnitudes_estimated_known, magnitudes_estimated_twostep, magnitudes_estimated_joint)
axs = (ax1, ax2, ax3)

for (ax, nbin, mag) in zip(axs, [20, 20, 40], mags)
    error = rm_na(magnitudes_true .- mag)
    hist!(ax, error, bins = nbin, strokecolor = :black, strokewidth = 1, label = "error")
    lines!(ax, repeat([0.0], 2), [0.0, 150], linestyle = :dash, color = :red, label = "zero error")
    #lines!(ax, repeat([Statistics.mean(error)], 2), [0.0, 150], linestyle = :dash, color = :orange, label = "mean error")
    #lines!(ax, repeat([Statistics.median(error)], 2), [0.0, 150], linestyle = :dash, color = :green, label = "median error")
end
linkxaxes!(axs...)
fig




fig = Figure(size = (600, 750))
ax1 = Axis(fig[1,1], xgridvisible = false, ygridvisible = false, title = L"\text{known }𝛌,𝛍,\text{ unknown } η", xscale = log10)
ax2 = Axis(fig[2,1], xgridvisible = false, ygridvisible = false, title = L"\text{unknown }𝛌,𝛍,η\text{, two step}", xscale = log10)
ax3 = Axis(fig[3,1], xgridvisible = false, ygridvisible = false, title = L"\text{unknown }𝛌,𝛍,η\text{, joint}", xscale = log10, 
        ylabel = L"\text{error} = \text{mag}_\text{true} - \text{mag}_\text{estimated}",
        xlabel = L"\text{number of tips}")

mags = (magnitudes_estimated_known, magnitudes_estimated_twostep, magnitudes_estimated_joint)
axs = (ax1, ax2, ax3)

for (ax, nbin, mag) in zip(axs, [20, 20, 40], mags)
    error = magnitudes_true .- mag
    scatter!(ax, ntip_all, error)
    #hist!(ax, error, bins = nbin, strokecolor = :black, strokewidth = 1, label = "error")
    lines!(ax, [extrema(ntip_all)...], repeat([0.0], 2), linestyle = :dash, color = :red, label = "zero error")
    #lines!(ax, repeat([Statistics.mean(error)], 2), [0.0, 150], linestyle = :dash, color = :orange, label = "mean error")
    #lines!(ax, repeat([Statistics.median(error)], 2), [0.0, 150], linestyle = :dash, color = :green, label = "median error")
end
linkxaxes!(axs...)
fig












magnitudes_true

magnitudes_estimated_joint


3+ 5













