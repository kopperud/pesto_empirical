using Revise
using Glob
using ProgressMeter
using CSV
using JLD2
using DataFrames
using RCall
using CairoMakie



datapaths_known = Glob.glob("output/simulations/up_vs_down/jld2/*.jld2")

datapaths_twostep = Glob.glob("output/simulations/up_vs_down_netdiv_unknown/jld2/*.jld2")

datapaths_joint = Glob.glob("output/simulations/up_vs_down_jointhyper/jld2/*.jld2")

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

N_estimated_known = zeros(Float64, ntrees, 5, 5)
for (i, path) in enumerate(datapaths_known)
    tree_index = parse(Int64, split(basename(path), ".")[1])
    N_estimated_known[tree_index,:,:] .= load(path, "Nsum")
end

Δrs_twostep = zeros(Float64, ntrees, 36, 36)
N_estimated_twostep = zeros(Float64, ntrees, 36, 36)
N_estimated_twostep[:,:,:] .= NaN
for (i, path) in enumerate(datapaths_twostep)
    tree_index = parse(Int64, split(basename(path), ".")[1])
    N_estimated_twostep[tree_index,:,:] .= load(path, "Nsum")

    λ = load(path, "λ")
    μ = load(path, "μ")
    r = λ .- μ
    Δr = r .- r'

    Δrs_twostep[tree_index,:,:] .= Δr
end

Δrs_joint = zeros(Float64, ntrees, 36, 36)
N_estimated_joint = zeros(Float64, ntrees, 36, 36)
N_estimated_joint[:,:,:] .= NaN
for (i, path) in enumerate(datapaths_joint)
    tree_index = parse(Int64, split(basename(path), ".")[1])
    N_estimated_joint[tree_index,:,:] .= load(path, "Nsum")

    λ = load(path, "λ")
    μ = load(path, "μ")
    r = λ .- μ
    Δr = r .- r'

    Δrs_joint[tree_index,:,:] .= Δr
end


size(N_estimated_joint)



#Δr
r = [0.01, 0.03, 0.05, 0.07, 0.09]
Δr = r .- r'
magnitudes_true = Float64[]
for i in 1:size(Ntrue)[1]
    mag = sum(Ntrue[i,:,:] .* Δr) / sum(Ntrue[i,:,:])
    if !isnan(mag)
        push!(magnitudes_true, mag)
    end
end

magnitudes_estimated_known = Float64[]
for i in 1:size(N_estimated_known)[1]
    mag = sum(N_estimated_known[i,:,:] .* Δr) / sum(N_estimated_known[i,:,:])
    if !isnan(mag)
        push!(magnitudes_estimated_known, mag)
    end
end


magnitudes_estimated_twostep = Float64[]
for i in 1:size(N_estimated_twostep)[1]
    Δr = Δrs_twostep[i,:,:]
    mag = sum(N_estimated_twostep[i,:,:] .* Δr) / sum(N_estimated_twostep[i,:,:])
    if !isnan(mag)
        push!(magnitudes_estimated_twostep, mag)
    end
end

magnitudes_estimated_joint = Float64[]
for i in 1:size(N_estimated_joint)[1]
    Δr = Δrs_joint[i,:,:]
    mag = sum(N_estimated_joint[i,:,:] .* Δr) / sum(N_estimated_joint[i,:,:])
    if !isnan(mag)
        push!(magnitudes_estimated_joint, mag)
    end
end



fig = Figure(size = (450, 750))

ax1 = Axis(fig[1,1], 
    title = L"\text{true shifts (reconstructed trees)}",
    ylabel = L"\text{number of trees}")
ax2 = Axis(fig[2,1], 
    title = L"\text{estimated shifts, known }𝛌,𝛍, \text{ unknown }η",
    ylabel = L"\text{number of trees}")
ax3 = Axis(fig[3,1], 
    title = L"\text{estimated shifts, unknown }𝛌,𝛍,η\text{ (two-step, }n^2 = 36)",
    xlabel = "magnitude",
    ylabel = L"\text{number of trees}")
ax4 = Axis(fig[4,1], 
        title = L"\text{estimated shifts, unknown }𝛌,𝛍,η\text{ (joint, }n^2 = 36)",
        xlabel = "magnitude",
        ylabel = L"\text{number of trees}")
linkaxes!(ax1, ax2, ax3)
linkxaxes!(ax1, ax4)
hist!(ax1, magnitudes_true, bins = 50)
hist!(ax2, magnitudes_estimated_known, bins = 50)
hist!(ax3, magnitudes_estimated_twostep, bins = 50)
hist!(ax4, magnitudes_estimated_joint, bins = 20)
fig

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

save("figures/Q2_simulations.pdf", fig)


## eta estimates
η1 = [load(path, "etaml") for path in datapaths_known]
η2 = [load(path, "etaml") for path in datapaths_twostep]
η3 = [load(path, "etaml") for path in datapaths_joint]

attempts = [load(path, "n_attempts") for path in datapaths_joint]


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




fig2 = Figure(size = (400, 800));

ax1 = Axis(fig2[1,1], xscale = log10, yscale = log10, title = "true", xticklabelsvisible = false, ylabel = "η", topspinevisible = false, rightspinevisible = false, xminorgridvisible = false, yminorgridvisible = false)
ax2 = Axis(fig2[2,1], xscale = log10, yscale = log10, title = "estimated (λ,μ known, η unkown)", xticklabelsvisible = false, ylabel = "η", topspinevisible = false, rightspinevisible = false, xminorgridvisible = false, yminorgridvisible = false)
ax3 = Axis(fig2[3,1], xscale = log10, yscale = log10, title = "estimated (λ,μ,η unknown), two-step", xticklabelsvisible = false, ylabel = "η", topspinevisible = false, rightspinevisible = false, xminorgridvisible = false, yminorgridvisible = false)
ax4 = Axis(fig2[4,1], xscale = log10, yscale = log10, title = "estimated (λ,μ,η unknown), joint", ylabel = "η", topspinevisible = false, rightspinevisible = false, xminorgridvisible = false, yminorgridvisible = false)

linkaxes!(ax1, ax2, ax3, ax4)

true_shift_rate = [sum(item) for item in items] ./ tree_lengths_all
true_shift_rate[true_shift_rate .== 0] .= 1e-8

scatter!(ax1, ntip_all, true_shift_rate, bins = 50)
scatter!(ax2, ntip1, η1, bins = 50)
scatter!(ax3, ntip2, η2, bins = 50)
scatter!(ax4, ntip3, η3, bins = 50)
for ax in (ax1, ax2, ax3, ax4)
    lines!(ax, [2.0, 10_000.0], [0.001, 0.001], color = :red, linestyle = :dash) 
end

fig2


