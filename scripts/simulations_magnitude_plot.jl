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
inference = "empirical"

df = CSV.read("output/age_scaling_effect_munged.csv", DataFrame)
#df = df[df[!,:inference] .== inference,:]

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
#names = unique(df[!,:name])
fpaths = Glob.glob("output/simulations/age_scaling_effect/jld2/*.jld2")
@showprogress for fpath in fpaths
    fname = split(Base.basename(fpath), ".")[1]
    height, i = split(fname, "_")
    
    height = parse(Float64, height[2:end])
    i = parse(Int64, i)

    x = JLD2.load(fpath)
    d[height, i] = x
end

#models = Dict{(Float64, Int64), SSEconstant}()
models = Dict()
heights = [30, 40, 50, 60, 70, 80, 90, 100]
iters = collect(1:500)
#for name in names
for height in heights
    for i in iters
        λ = d[height, i]["lambda"]
        μ = d[height, i]["mu"]
        η = d[height, i]["etaml"]
        models[height, i] = SSEconstant(λ, μ, η)
    end
end


#trees = Dict()
datasets = Dict()
@showprogress for height in heights
    for i in iters
        fpath = string("data/simulations/age_scaling_effect/h", height, "_", i, ".tre")
        tree = readtree(fpath)
        
        if all(tree.edge_length .> 0)
            #bn = Base.Filesystem.basename(fpath)
            #trees[bn] = tree
            datasets[height,i] = SSEdata(tree, 1.0) ## complete taxon sampling
        end
    end
end


rates = Dict()
for height in heights
    for i in iters
        fpath = string("output/simulations/age_scaling_effect/rates/h", height, "_", i, ".csv")
        rates_df = CSV.read(fpath, DataFrame)
        sort!(rates_df, :edge)
        rates[height, i] = rates_df[2:end,:]
    end
end

is_supported = Dict()
for height in heights
    for i in iters
        rdf = rates[height, i]
        iss = ((rdf[!,:nshift] .> 0.5) .& (rdf[!,:shift_bf] .> 10))
        is_supported[height, i] = iss
    end
end

##################
##
##     calculate magnitudes
##
###############

#n_datasets = length(datasets)
#n_datasets = length(d)
#heights = zeros(n_datasets)
magnitudes = zeros(length(heights), length(iters), 2)


#for (i, name) in enumerate(keys(d))
#heights_vector = 
for (height_index, height) in enumerate(heights)
    for i in iters
        model = models[height,i]
        #heights[i] = maximum(datasets[name .* ".tree"].node_depth)
        
        Nsum = sum(d[height,i]["N"], dims = 1)[1,:,:]
        m = magnitude(model, Nsum)
        magnitudes[height_index, i, 1] = m

        support = is_supported[height,i]
        Nsum = sum(d[height,i]["N"][support,:,:], dims = 1)[1,:,:]
        m = magnitude(model, Nsum)
        magnitudes[height_index, i, 2] = m
    end 
end

plotdf = DataFrame(
    "magnitudes_pooled" => vcat(magnitudes[:,:,1]...),
    "magnitudes_supported" => vcat(magnitudes[:,:,2]...),
    "heights" => repeat(heights, 500),
    "replicate" => repeat(collect(1:500), 8)
    #"name" => collect(keys(d))
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


fig2 = Figure(size = (650, 300), fontsize = 14, 
                figure_padding = (5,8,1,1))

magnitudes_pooled = plotdf[!,:magnitudes_pooled]
magnitudes_support = support_df[!, :magnitudes_supported]
heights_pooled = plotdf[!,:heights]
heights_support = support_df[!,:heights]

xt = collect(range(extrema(plotdf[!,:magnitudes_pooled])...; length = 5))
xtl = [@sprintf("%.2f", x) for x in xt]

ax_hist_pooled = Axis(fig2[1,1], 
            ylabel = L"\text{frequency}", 
            xlabel = L"\text{magnitude }(\Delta r)",
            title = L"\text{all branches pooled}",
            xgridvisible = false, 
            ygridvisible = false,
            xticks = (xt, xtl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelsize = 9,
            yticklabelsize = 9)

xt = collect(range(extrema(support_df[!,:magnitudes_supported])...; length = 5))
xtl = [@sprintf("%.2f", x) for x in xt]

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


CairoMakie.hist!(ax_hist_pooled, 
    support_df[!,:magnitudes_pooled], bins = 10, color = "gray")
CairoMakie.hist!(ax_hist_supported, 
    support_df[!,:magnitudes_supported], bins = 10, color = "gray")


xt = 10 .^ (collect(range(extrema(log10.(plotdf[!,:heights]))...; length = 5)))
xtl = [@sprintf("%.1f", x) for x in xt]


if mag_logscale
    yr = collect(range(extrema(log10.(plotdf[!,:magnitudes_pooled]))...; length = 5))
    yt = 10 .^ (yr)
else
    yr = collect(range(extrema(plotdf[!,:magnitudes_pooled])...; length = 5))
    yt = yr
end
ytl = [@sprintf("%.2f", y) for y in yt]

ax_scatter_pooled = Axis(fig2[2,1], 
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
    yt = 10 .^ (yr)
    yr = collect(range(extrema(log10.(support_df[!,:magnitudes_supported]))...; length = 5))
else
    yr = collect(range(extrema(support_df[!,:magnitudes_supported])...; length = 5))
    yt = yr
end
ytl = [@sprintf("%.2f", y) for y in yt]

ax_scatter_support = Axis(fig2[2,2], 
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



#for (_df, ax) in zip([plotdf, support_df], [ax_scatter_pooled, ax_scatter_support])
for (magnitudes, heights, ax) in zip([magnitudes_pooled, magnitudes_support], [heights_pooled, heights_support], [ax_scatter_pooled, ax_scatter_support])

    if mag_logscale
        β, Varβ, ySE = ols_regression(log10.(heights), log10.(magnitudes))
    else
        β, Varβ, ySE = ols_regression(log10.(heights), magnitudes)
    end

    x = collect(Pesto.lrange(30.0, 100.0, 20))
    linefit(x) = β[1] + β[2]*log10.(x)

    if mag_logscale
        ylog = linefit.(x)
    else
        yarithmetic = linefit.(x)
    end

    yVar = Varβ[1,1] .+ 2 .* log10.(x) .* Varβ[1,2] .+ (log10.(x) .^ 2) .* Varβ[2,2]

    if mag_logscale
        yupper = 10 .^(ylog .+ 2*sqrt.(yVar))
        ylower = 10 .^(ylog .- 2*sqrt.(yVar))
        y = 10 .^ linefit.(x)
    else
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

colgap!(fig2.layout, 5)
rowgap!(fig2.layout, 5)
fig2

CairoMakie.save("figures/magnitude_simulatedtrees_empiricalbayes.pdf", fig2)


for (i, name) in enumerate(keys(d))
    println(name, ": \t", magnitudes[i])
end


ntips = Int64[]
heights_vector = Float64[]
mags = Float64[]

for (height_index, height) in enumerate(heights)
    for i in iters
        tn = length(datasets[height, i].tiplab)
        append!(ntips, tn)

        h = float(height)
        append!(heights_vector, h)

        append!(mags, magnitudes[height_index,i,1])
    end
end



#fig3 = Figure()
#xt = 
ax = Axis(fig2[1:2,3], 
    xlabel = L"\text{number of tips}", 
    ylabel = L"\text{magnitude }(\Delta r)",
    xgridvisible = false, 
    ygridvisible = false,
    xscale = log10,
    title = L"\text{all branched pooled}",
    topspinevisible = false,
    rightspinevisible = false)
CairoMakie.scatter!(ax, ntips, mags, color = "black")
CairoMakie.lines!(ax, [2, 30_000], [0.0, 0.0], color = "red", linestyle = :dash)


fig2

CairoMakie.save("figures/magnitude_simulatedtrees_empiricalbayes.pdf", fig2)


fig3
ntips



