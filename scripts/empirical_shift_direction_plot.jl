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

function direction(model, N)
    r = model.λ .- model.μ

    Δr = r .- r'

    Δr_sign = sign.(Δr)

    mean_direction = sum(Δr_sign .* N) / sum(N)
    return(mean_direction)
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

##################
##
##     calculate direction
##
###############

#n_datasets = length(datasets)
n_datasets = length(d)
directions = zeros(n_datasets)
heights = zeros(n_datasets)

for (i, name) in enumerate(keys(d))
    model = models[name]
    heights[i] = maximum(datasets[name .* ".tree"].node_depth)
    
    Nsum = sum(d[name]["N"], dims = 1)[1,:,:]
    dir = direction(model, Nsum)

    directions[i] = dir
end

plotdf = DataFrame(
    "directions" => directions,
    "heights" => heights,
    "name" => collect(keys(d))
)

##################
##
##   set up the makie figure
##
###############

# fig1 = Figure(size = (450, 180), fontsize = 14, figure_padding = (5,8,1,1));
fig1 = Figure(size = (220, 180), fontsize = 14, figure_padding = (5,8,1,1));

xt = collect(range(extrema(directions)...; length = 5))
xt = [-0.60, -0.30, 0.0, 0.30, 0.60, 0.90, 1.2]
xtl = [@sprintf("%.2f", x) for x in xt]

ax_hist = Axis(fig1[1,1], 
            ylabel = L"\text{frequency}", 
            xlabel = L"\text{shift direction}",
            xgridvisible = false, 
            ygridvisible = false,
            xticks = (xt, xtl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelsize = 9,
            yticklabelsize = 9)

CairoMakie.hist!(ax_hist, 
    directions, bins = 14, color = "gray")

#CairoMakie.save("figures/directions_empirical_histogram.pdf", fig1)

xt = 10 .^ (collect(range(extrema(log10.(plotdf[!,:heights]))...; length = 5)))
xtl = [@sprintf("%.1f", x) for x in xt]


yr = collect(range(extrema(plotdf[!,:directions])...; length = 5))
yr = [-0.60, -0.30, 0.0, 0.30, 0.60, 0.90, 1.2]
yt = yr
ytl = [@sprintf("%.2f", y) for y in yt]


ax_scatter = Axis(fig1[1,2], 
        ylabel = L"\text{shift direction}", 
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



yr = collect(range(extrema(support_df[!,:directions])...; length = 5))
yt = yr
ytl = [@sprintf("%.2f", y) for y in yt]

β, Varβ, ySE = ols_regression(log10.(heights), directions)

x = collect(Pesto.lrange(extrema(heights)..., 20))
linefit(x) = β[1] + β[2]*log10.(x)
yarithmetic = linefit.(x)

#yVar = Varβ[1,1] .+ 2 .* x .* Varβ[1,2] .+ (x .^ 2) .* Varβ[2,2]
yVar = Varβ[1,1] .+ 2 .* log10.(x) .* Varβ[1,2] .+ (log10.(x) .^ 2) .* Varβ[2,2]
yupper = yarithmetic .+ 2*sqrt.(yVar)
ylower = yarithmetic .- 2*sqrt.(yVar)
y = linefit.(x)

CairoMakie.band!(ax_scatter, x, ylower, yupper, color = "#e0e0e0")
CairoMakie.lines!(ax_scatter, x, y; label = "OLS", markersize = 7, color = "gray", linestyle = :dash)

CairoMakie.scatter!(ax_scatter, 
                    heights, 
                    directions,
                    color = "black",
                    markersize = 7)

#lines!(ax_scatter_pooled, [extrema(heights)...], [0.0, 0.0], linestyle = :dash, color = :red)
colgap!(fig1.layout, 5)
rowgap!(fig1.layout, 5)
fig1



###################
n_positive = sum(directions .> 0)
n_negative = sum(directions .< 0)

println("$n_positive directions are positive, while $n_negative are negative")

#CairoMakie.save("figures/directions_empirical_joint.pdf", fig1)

for (i, name) in enumerate(keys(d))
    println(name, ": \t", directions[i])
end
