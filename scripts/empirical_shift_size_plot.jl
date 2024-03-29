using Distributions
using Glob, DataFrames, CSV, RCall, ProgressMeter
using LaTeXStrings, Measures
using JLD2
using Printf
using CairoMakie
using Pesto

########################
##
##   some helper functions
##
##########################

function ols_regression(x, y)
    X = hcat([1 for _ in 1:length(x)], x)
    n, p = size(X)

    ## OLS
    β = (X' * X) \ X' * y
    yhat = X * β
    sigma_squared = (1 / (n - p - 1)) * sum((y .- yhat).^2) ## MLE for sigma^2
    #s_squared = (y .- yhat)' * (y .- yhat) ./ (n - p) ## OLS for sigma^2
    Varβ = inv(X' * X) .* sigma_squared
    yVar = x -> Varβ[1,1] + 2*x*Varβ[1,2] + (x^2)*Varβ[2,2]
    ySE = x -> sqrt.(yVar(x))
    return(β, Varβ, ySE)
end

function support_vector(name, subdir)
    df = CSV.read("output/" * subdir * "/rates/" * name * ".csv", DataFrame)
    sort!(df, [:edge])
    df = df[2:end,:]

    is_supported = Bool[]
    for row in eachrow(df)
        if (row[:nshift] > 0.5) & (row[:shift_bf] > 10)
            s = true
        else
            s = false
        end
        append!(is_supported, s)
    end

    return(is_supported)
end

function compute_ratios(name, model, subdir, filter = "")
    N = d[name]["N"]

    if filter == "support"
        is_supported = support_vector(name, subdir)
        N = N[is_supported,:,:]
    end

    Nmatrix = sum(N, dims = 1)[1,:,:]
    nbins = 14

    Ns = zeros(3,nbins)
    filters = ["extinction", "speciation", ""]
    limits = [-1.2, 1.2]

    dfs = []
    for (i, filter) in enumerate(filters)
        mids, bins = makebins(Nmatrix, model, limits...; filter = filter, nbins = nbins)
        Ns[i,:] = bins[:,3]
        df = DataFrame(
            "Δr" => bins[:,3],
            "mids" => mids,
            "subset" => [i for _ in 1:nbins]
        )
        append!(dfs, [df])
    end
    df = DataFrame(
        "Δr" => dfs[3][!, "Δr"] .- dfs[1][!, "Δr"] .- dfs[2][!, "Δr"],
        "mids" => dfs[3][!, "mids"],
        "subset" => 3
    )
    dfs[3] = df

    Nλ = sum(dfs[1][!,:Δr])
    Nμ = sum(dfs[2][!,:Δr])
    Njoint = sum(dfs[3][!,:Δr])
    Nall = sum([Nλ, Nμ, Njoint])
    rs = [Nλ, Nμ, Njoint] ./ Nall
    return(rs)
end


#inference = "empirical_fixedprior"
inference = "empirical"

df = CSV.read("output/empirical_munged.csv", DataFrame)
df = df[df[!,:inference] .== inference,:]

## N tensor
#
# Dimensions:
# 1: edge index (M)
# 2: row in the N matrix (index i), arrival state
# 3: column in the N matrix (index j), departure state
d = Dict()
names = unique(df[!,:name])
fpaths = Glob.glob("output/" * inference * "/jld2/*.jld2")
#fpaths = Glob.glob("output/empirical/jld2/*.jld2")
@showprogress for fpath in fpaths
    name = split(Base.basename(fpath), ".")[1]
    x = JLD2.load(fpath)
    d[name] = x
end

models = Dict{String, SSEconstant}()
for name in names
    λ = d[name]["lambda"]
    μ = d[name]["mu"]
    η = d[name]["etaml"]
    models[name] = SSEconstant(λ, μ, η)
end



#fpaths = glob("*.tree", "data/empirical")
fpaths = ["data/empirical/" * name * ".tree" for name in names]

df = CSV.read("data/empirical/Phylogenies for DeepILS project.csv", DataFrame)
ρs = Dict()
for row in eachrow(df)
    fn = row["Filename"]
    ρ = row["P extant sampling"]
    ρs[fn] = ρ
end

trees = Dict()
datasets = Dict()
for fpath in fpaths
    println(fpath)
    tree = readtree(fpath)
    
    if all(tree.edge_length .> 0)
        bn = Base.Filesystem.basename(fpath)
        trees[bn] = tree
        datasets[bn] = SSEdata(trees[bn], ρs[bn])
    end
end

#heights = [maximum(d.node_depth) for (key, d) in datasets]




####################################
##
##   Fig. 1: The rate shift type (speciation vs extinction)
##
####################################
fig = Figure(size = (550, 220), fontsize = 14, 
            figure_padding = (1,1,1,1))

colors = [:steelblue, "orange", "gray"]
labels = [
    L"\text{Shift in \lambda}", 
    L"\text{Shift in \mu}", 
    L"\text{Shift in both}"
    ]


g = fig[1,1] = GridLayout()

name_subset = [
    "Primates_Springer2012",
    "Mammalia_AlvarezCarretero2022",
    "Rosidae_Sun2020",
    "Actinopterygii_Rabosky2018",
    "Agaricomycetes_SanchezGarcia2020"
]

titles = []
dnames = collect(keys(datasets))
for dname in dnames
    append!(titles, [split(dname, "_")[1]])
end


#datasets_order = [34, 9, 28, 29, 4]

axs = []
q = 1
for i in 1:2, j in 1:3
    if (i == 1) & (j == 3)
        #ax = Axis(fig[i,j], scene = false)
        #hidedecorations!(ax)
    else
        if i < 2
            xlabel = ""
        else
            xlabel = L"\Delta r"
        end

        title = split(name_subset[q], "_")[1]

        #foo_idx = datasets_order[q]
        ax = Axis(g[i,j], 
        xgridvisible = false, 
        ygridvisible = false,
        #title = titles[foo_idx],
        title = title,
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9)

        if i < 2
            hidexdecorations!(ax, ticks = false)
        end
    
        if j > 1
            hideydecorations!(ax, ticks = false, ticklabels = false)
        end

        append!(axs, [ax])
        q += 1
    end
end



netdiv_extrema = [-1.2, 1.2]

ratios = zeros(length(name_subset), 3)

for (q, name) in enumerate(name_subset)
    model = models[name]

    #Nmatrix = N[dataset_index,:,:,1]
    Nmatrix = sum(d[name]["N"], dims = 1)[1,:,:]
    nbins = 14

    Ns = zeros(3,nbins)
    filters = ["extinction", "speciation", ""]

    dfs = []
    for (i, filter) in enumerate(filters)
        mids, bins = makebins(Nmatrix, model, netdiv_extrema...; filter = filter, nbins = nbins)
        Ns[i,:] = bins[:,3]
        df = DataFrame(
            "Δr" => bins[:,3],
            "mids" => mids,
            "subset" => [i for _ in 1:nbins]
        )
        append!(dfs, [df])
    end
    df = DataFrame(
        "Δr" => dfs[3][!, "Δr"] .- dfs[1][!, "Δr"] .- dfs[2][!, "Δr"],
        "mids" => dfs[3][!, "mids"],
        "subset" => 3
    )
    dfs[3] = df
    
    Nλ = sum(dfs[1][!,:Δr])
    Nμ = sum(dfs[2][!,:Δr])
    Njoint = sum(dfs[3][!,:Δr])
    Nall = sum([Nλ, Nμ, Njoint])
    rs = [Nλ, Nμ, Njoint] ./ Nall
    ratios[q,:] .= rs



    bar_df = vcat(dfs...)

    barplot!(axs[q], 
        bar_df[!, "mids"], 
        bar_df[!, "Δr"], 
        stack = bar_df[!, "subset"],
        color = [colors[x] for x in bar_df[!, "subset"]],
        )
end

elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
title = "Shift in"
        #glegend = 
leg = Legend(g[1,3], elements, labels, labelsize=9, 
            patchsize = (10.0f0, 10.0f0),
            framevisible = false)




leg.tellheight = true
leg.tellwidth = true

fig
#ylabel = Label(g[1:2, 1], "Now the label is centered!", rotation = π/2)
fig

linkxaxes!(axs...)
ylabel = Label(g[1:2, 0], L"\text{number of rate shifts }(\hat{N})", rotation = π/2)
xlabel = Label(g[3, 1:3], L"\text{shift size in net diversification }(\Delta r)")


##############
#
#   VIOLIN FIGURE 
#
############


ratios1 = zeros(length(datasets), 3)
ratios2 = zeros(length(datasets), 3)

for (dataset_index, name) in enumerate(names)
    model = models[name]

    #dname = replace(name, "Toninietal2016" => "Toninietal.2016")
    #dname = replace(name, "Zanneetal2014rooteddated" => "Zanne.et.al.2014rooted.dated")


    rs1 = compute_ratios(name, 
            models[name],
            inference,
            "")

    rs2 = compute_ratios(name, 
            models[name],
            inference,
            "support")


    ratios1[dataset_index,:] .= rs1
    ratios2[dataset_index,:] .= rs2
end

## filter ratios for NaN
ratios2 = transpose(hcat([row for row in eachrow(ratios2) if !any(isnan.(row))]...))


ax_scatter1 = Axis(g[1:2,4], 
        ylabel = L"N_\text{rate} / N_\text{all}", 
        xgridvisible = false, 
        ygridvisible = false,
        xticks = (1:3, [L"\lambda", L"\mu", L"\text{both}"]),
        #title = L"\text{all branches}",
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelsize = 9,
        yticklabelsize = 9)

#= ax_scatter2 = Axis(g[1:2,5], 
        ylabel = "", 
        xgridvisible = false, 
        ygridvisible = false,
        xticks = (1:3, [L"\lambda", L"\mu", L"\text{both}"]),
        title = L"\text{supported branches}",
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelsize = 9,
        yticklabelsize = 9) =#

#= hideydecorations!(ax_scatter2, ticks = false) =#

## hard code the priors
priors = [1/11, 1/11, 9/11]
prior_ratios = zeros(length(datasets), 3)
for i in 1:size(prior_ratios)[1]
    prior_ratios[i,:] .= priors
end

is_larger1 = ratios1 .> prior_ratios
#is_larger2 = ratios1 .> prior_ratios
sum(is_larger1, dims = 1)
#sum(is_larger2, dims = 1)

for (i, prior) in enumerate(priors)
    x = [i-0.5, i+0.5]
    y = [prior, prior]
    CairoMakie.lines!(ax_scatter1, x, y, color = colors[i], linestyle = :solid, alpha = 0.5, label = "prior")
    #CairoMakie.lines!(ax_scatter2, x, y, color = colors[i], linestyle = :solid, alpha = 0.5, label = "prior")
end

xs = vcat([repeat([i], size(ratios1)[1]) for i in 1:3]...)
cs = vcat([repeat([i], size(ratios1)[1]) for i in colors]...)

CairoMakie.rainclouds!(ax_scatter1, xs .- 0.25, vcat(ratios1...), 
                     color = cs, label = "Posterior",
                     clouds=nothing, boxplot_width=0.5,
                     jitter_width=0.35,
                     side_nudge = 0.5,
                     markersize=4)

#= xs = vcat([repeat([i], size(ratios2)[1]) for i in 1:3]...)
cs = vcat([repeat([i], size(ratios2)[1]) for i in colors]...)

CairoMakie.rainclouds!(ax_scatter2, xs .- 0.25, vcat(ratios2...), 
                     color = cs, label = "Posterior",
                     clouds=nothing, boxplot_width=0.5,
                     jitter_width=0.35,
                     side_nudge = 0.5,
                     markersize=4) =#

xlabel2 = Label(g[3, 4], L"\text{posterior vs prior}")
#axislegend(ax_scatter1)
fig

for col in 1:3
    colsize!(g, col, Relative(0.20))
end
colsize!(g, 0, Relative(0.04))
colsize!(g, 4, Relative(0.35))
#colsize!(g, 5, Relative(0.22))

for i in 1:2
    rowsize!(g, i, Relative(0.45))
end
rowsize!(g, 3, Relative(0.1))

colgap!(g, 5)
rowgap!(g, 3)


#linkaxes!(ax_scatter1, ax_scatter2)

fig
#set_theme!(fig, figure_padding = 0)
#CairoMakie.save("figures/fig1_fixedprior.pdf", fig)
CairoMakie.save("figures/fig1_empiricalbayes.pdf", fig)
#CairoMakie.save("figures/histogram_magnitude_etaoneshift.pdf", fig)


