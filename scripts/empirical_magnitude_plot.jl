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

function magnitude(model, N)
    r = model.λ .- model.μ
    Δr = N .* (r .- r')
    mean_magnitude = sum(Δr) / sum(N)
    return(mean_magnitude)
end

function var_bds(model, N)
    r = (model.λ .- model.λ') .- (model.μ .- model.μ')

    m = magnitude(model, N)
    K = size(model.λ)[1]

    res = 0.0
    for i in 1:K, j in 1:K
        res += (r[i,j] * N[i,j] - m)^2
    end
    res = res / sum(N)
    return(res)
end

########################
##
##         read files
##
##########################

#inference = "empirical_fixedprior"
inference = "empirical_joint"

df = CSV.read("output/empirical_munged.csv", DataFrame)
df = df[df[!,:inference] .== inference,:]

########################
##
##           N tensor
## Dimensions:
## 1: edge index (M)
## 2: row in the N matrix (index i), arrival state
## 3: column in the N matrix (index j), departure state
##
##########################


d = Dict()
names = unique(df[!,:name])
fpaths = Glob.glob("output/" * inference * "/jld2/*.jld2")
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

fpaths = Glob.glob("output/" * inference * "/rates/*.csv")
rates = Dict{String, DataFrame}()
for fpath in fpaths
    name = split(Base.basename(fpath), ".")[1]
    rates_df = CSV.read(fpath, DataFrame)
    sort!(rates_df, :edge)
    rates[name] = rates_df[2:end,:]
end

is_supported = Dict{String, BitVector}()
for (name, rdf) in rates
    iss = ((rdf[!,:nshift] .> 0.5) .& (rdf[!,:shift_bf] .> 10))
    is_supported[name] = iss
end

##################
##
##     calculate magnitudes
##
###############

#n_datasets = length(datasets)
n_datasets = length(d)
magnitudes = zeros(n_datasets,2)
heights = zeros(n_datasets)

for (i, name) in enumerate(keys(d))
    model = models[name]
    heights[i] = maximum(datasets[name .* ".tree"].node_depth)
    
    support = is_supported[name]
    Nsum = sum(d[name]["N"][support,:,:], dims = 1)[1,:,:]
    m = magnitude(model, Nsum)

    magnitudes[i,2] = m

    Nsum = sum(d[name]["N"], dims = 1)[1,:,:]
    m = magnitude(model, Nsum)
    magnitudes[i,1] = m    
end

plotdf = DataFrame(
    "magnitudes_pooled" => magnitudes[:,1],
    "magnitudes_supported" => magnitudes[:,2],
    "heights" => heights,
    "name" => collect(keys(d))
)
support_df = deepcopy(plotdf)
filter!(:magnitudes_supported => x -> !isnan(x), support_df)

#magnitudes = plotdf[!,:magnitudes]
#heights = plotdf[!,:heights]
#n_datasets = size(plotdf)[1]

## some variance bs, don't remember why I did this
#= variances = zeros(n_datasets)
for (i, name) in enumerate(keys(d))
    model = models[name]
    Nsum = sum(d[name]["N"], dims = 1)[1,:,:]

    variances[i] = var_bds(model, Nsum)
end
CairoMakie.scatter(log10.(heights), log10.(variances)) =#

##################
##
##   set up the makie figure
##
###############
mag_logscale = false


fig2 = Figure(size = (450, 180), fontsize = 14, 
                figure_padding = (5,8,1,1));

magnitudes_pooled = plotdf[!,:magnitudes_pooled]
magnitudes_support = support_df[!, :magnitudes_supported]
heights_pooled = plotdf[!,:heights]
heights_support = support_df[!,:heights]

#xt = collect(range(extrema(magnitudes_pooled)...; length = 5))
xt = [-0.25, 0.0, 0.25, 0.50, 0.75, 1.0, 1.25]
xtl = [@sprintf("%.2f", x) for x in xt]

ax_hist_pooled = Axis(fig2[1,1], 
            ylabel = L"\text{frequency}", 
            xlabel = L"\text{magnitude }(\Delta r)",
            #title = L"\text{all branches pooled}",
            xgridvisible = false, 
            ygridvisible = false,
            xticks = (xt, xtl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelsize = 9,
            yticklabelsize = 9)

#xt = collect(range(extrema(support_df[!,:magnitudes_supported])...; length = 5))
#xt = [-0.25, 0.0, 0.25, 0.50, 0.75, 1.0, 1.25]
#xtl = [@sprintf("%.2f", x) for x in xt]

#= 
ax_hist_supported = Axis(fig2[1,2], 
            #ylabel = L"\text{frequency}", 
            xlabel = L"\text{magnitude }(\Delta r)",
            title = L"\text{supported branches}",
            xgridvisible = false, 
            ygridvisible = false,
            xticks = (xt, xtl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelsize = 9,
            yticklabelsize = 9)
CairoMakie.hist!(ax_hist_supported, 
    support_df[!,:magnitudes_supported], bins = 10, color = "gray")
 =#

CairoMakie.hist!(ax_hist_pooled, 
    magnitudes_pooled, bins = 15, color = "gray")


xt = 10 .^ (collect(range(extrema(log10.(plotdf[!,:heights]))...; length = 5)))
xtl = [@sprintf("%.1f", x) for x in xt]


if mag_logscale
    yr = collect(range(extrema(log10.(plotdf[!,:magnitudes_pooled]))...; length = 5))
    yt = 10 .^ (yr)
else
    yr = collect(range(extrema(plotdf[!,:magnitudes_pooled])...; length = 5))
    yr = [-0.25, 0.0, 0.25, 0.50, 0.75, 1.0, 1.25]
    yt = yr
end
ytl = [@sprintf("%.2f", y) for y in yt]

ax_scatter_pooled = Axis(fig2[1,2], 
        ylabel = L"\text{magnitude } (\Delta r)", 
        xgridvisible = false, 
        ygridvisible = false,
        #yscale = log10, 
        xscale = log10,
        xticks = (xt, xtl),
        yticks = (yt, ytl),
        xlabel = L"\text{tree height (Ma)}",
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelsize = 9,
        yticklabelsize = 9)


if mag_logscale
    yr = collect(range(extrema(log10.(support_df[!,:magnitudes_supported]))...; length = 5))
    yt = 10 .^ (yr)
else
    yr = collect(range(extrema(support_df[!,:magnitudes_supported])...; length = 5))
    yt = yr
end
ytl = [@sprintf("%.2f", y) for y in yt]

#= ax_scatter_support = Axis(fig2[2,2], 
        #ylabel = L"\text{magnitude } (\Delta r)", 
        xgridvisible = false, 
        ygridvisible = false,
        #yscale = log10, 
        xscale = log10,
        xticks = (xt, xtl),
        yticks = (yt, ytl),
        xlabel = L"\text{tree height (Ma)}",
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelsize = 9,
        yticklabelsize = 9)
 =#


#for (_df, ax) in zip([plotdf, support_df], [ax_scatter_pooled, ax_scatter_support])
#for (magnitudes, heights, ax) in zip([magnitudes_pooled, magnitudes_support], [heights_pooled, heights_support], [ax_scatter_pooled, ax_scatter_support])
for (magnitudes, heights, ax) in zip([magnitudes_pooled], [heights_pooled], [ax_scatter_pooled])
    if mag_logscale
        β, Varβ, ySE = ols_regression(log10.(heights), log10.(magnitudes))
    else
        β, Varβ, ySE = ols_regression(log10.(heights), magnitudes)
    end

    x = collect(Pesto.lrange(extrema(heights)..., 20))
    linefit(x) = β[1] + β[2]*log10.(x)

    if mag_logscale
        ylog = linefit.(x)
    else
        yarithmetic = linefit.(x)
    end

    

    if mag_logscale
        yVar = Varβ[1,1] .+ 2 .* log10.(x) .* Varβ[1,2] .+ (log10.(x) .^ 2) .* Varβ[2,2]
        yupper = 10 .^(ylog .+ 2*sqrt.(yVar))
        ylower = 10 .^(ylog .- 2*sqrt.(yVar))
        y = 10 .^ linefit.(x)
    else
        #yVar = Varβ[1,1] .+ 2 .* x .* Varβ[1,2] .+ (x .^ 2) .* Varβ[2,2]
        yVar = Varβ[1,1] .+ 2 .* log10.(x) .* Varβ[1,2] .+ (log10.(x) .^ 2) .* Varβ[2,2]
        yupper = yarithmetic .+ 2*sqrt.(yVar)
        ylower = yarithmetic .- 2*sqrt.(yVar)
        y = linefit.(x)
    end

    CairoMakie.band!(ax, x, ylower, yupper, color = "#e0e0e0")
    CairoMakie.lines!(ax, x, y; label = "OLS", markersize = 7, color = "gray", linestyle = :dash)

    CairoMakie.scatter!(ax, 
                        heights, 
                        magnitudes,
                        color = "black",
                        markersize = 7)
end

#lines!(ax_scatter_pooled, [extrema(heights)...], [0.0, 0.0], linestyle = :dash, color = :red)

colgap!(fig2.layout, 5)
rowgap!(fig2.layout, 5)
fig2

###################

n_positive = sum(magnitudes_pooled .> 0)
n_negative = sum(magnitudes_pooled .< 0)

println("$n_positive magnitudes are positive, while $n_negative are negative")


#CairoMakie.save("figures/magnitude_empiricalbayes_log.pdf", fig2)
#CairoMakie.save("figures/magnitude_empiricalbayes.pdf", fig2)
#CairoMakie.save("figures/magnitude_fixedprior.pdf", fig2)
#CairoMakie.save("figures/magnitude_empirical_joint.pdf", fig2)

for (i, name) in enumerate(keys(d))
    println(name, ": \t", magnitudes[i])
end

dfm = DataFrame(
    "name" => collect(keys(d)),
    "m" => magnitudes
)

sort(dfm, :m)


alphabet = "abcdefgh"
alphabet = "ABCDEFGH"
alphabet = ["i", "ii", "iii", "iv", "v"]

for i in 1:5
    for j in 1:5
        arrival   = alphabet[i]
        departure = alphabet[j]
        print("(\\lambda_{$arrival}, \\mu_{$departure})  ")

        if j < 5
            print(" & ")
        end
    end
    print("\\\\ \n")
end