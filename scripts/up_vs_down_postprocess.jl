using Revise
using Glob
using ProgressMeter
using CSV
using JLD2
using DataFrames
using RCall
using CairoMakie



datapaths = Glob.glob("output/simulations/up_vs_down/jld2/*.jld2")

datapaths2 = Glob.glob("output/simulations/up_vs_down_netdiv_unknown/jld2/*.jld2")

treepaths = Glob.glob("data/simulations/up_vs_down/*.tre")



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

ntrees = length(treepaths)
Ntrue = zeros(Int64, ntrees, 5, 5)
for tree_index in 1:ntrees
    Ntrue[tree_index,:,:] .= items[tree_index]
end

N_estimated = zeros(Float64, ntrees, 5, 5)
for i in eachindex(datapaths)
    tree_index = parse(Int64, split(basename(datapaths[i]), ".")[1])
    N_estimated[tree_index,:,:] .= load(datapaths[i], "Nsum")
end

Δrs = zeros(Float64, ntrees, 36, 36)
N_estimated_netdiv_unknown = zeros(Float64, ntrees, 36, 36)
N_estimated_netdiv_unknown[:,:,:] .= NaN
for i in eachindex(datapaths2)
    tree_index = parse(Int64, split(basename(datapaths2[i]), ".")[1])
    N_estimated_netdiv_unknown[tree_index,:,:] .= load(datapaths2[i], "Nsum")

    λ = load(datapaths2[i], "λ")
    μ = load(datapaths2[i], "μ")
    r = λ .- μ
    Δr = r .- r'

    Δrs[tree_index,:,:] .= Δr
end

size(N_estimated_netdiv_unknown)



Δr
r = [0.01, 0.03, 0.05, 0.07, 0.09]
Δr = r .- r'
#magnitudes_true = zeros(ntrees)
magnitudes_true = Float64[]
for i in 1:size(Ntrue)[1]
    mag = sum(Ntrue[i,:,:] .* Δr) / sum(Ntrue[i,:,:])
    if !isnan(mag)
        push!(magnitudes_true, mag)
    end
end

magnitudes_estimated = Float64[]
for i in 1:size(N_estimated)[1]
    mag = sum(N_estimated[i,:,:] .* Δr) / sum(N_estimated[i,:,:])
    if !isnan(mag)
        push!(magnitudes_estimated, mag)
    end
end


magnitudes_estimated_unknown = Float64[]
for i in 1:size(N_estimated)[1]
    Δr = Δrs[i,:,:]
    mag = sum(N_estimated_netdiv_unknown[i,:,:] .* Δr) / sum(N_estimated_netdiv_unknown[i,:,:])
    if !isnan(mag)
        push!(magnitudes_estimated_unknown, mag)
    end
end


fig = Figure(size = (450, 600))

ax1 = Axis(fig[1,1], 
    title = L"\text{true shifts (reconstructed trees)}",
    ylabel = L"\text{number of trees}")
ax2 = Axis(fig[2,1], 
    title = L"\text{estimated shifts, known }𝛌,𝛍, \text{ unknown }η",
    ylabel = L"\text{number of trees}")
ax3 = Axis(fig[3,1], 
    title = L"\text{estimated shifts, unknown }𝛌,𝛍,η\text{ (empirical Bayes, }n^2 = 36)",
    xlabel = "magnitude",
    ylabel = L"\text{number of trees}")
linkaxes!(ax1, ax2, ax3)
hist!(ax1, magnitudes_true, bins = 50)
hist!(ax2, magnitudes_estimated, bins = 50)
hist!(ax3, magnitudes_estimated_unknown, bins = 50)
fig

save("figures/Q2_simulations.pdf", fig)


