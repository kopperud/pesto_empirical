using Distributions
using Glob, DataFrames, CSV, RCall, ProgressMeter
using LaTeXStrings, Measures
using JLD2
using Printf
using CairoMakie
using Pesto

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


#inference = "empirical_fixedprior"
inference = "empirical"

df = CSV.read("output/empirical_munged.csv", DataFrame)
df = df[df[!,:inference] .== inference,:]
#df = df[df[!,:inference] .== "empirical",:]

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


#d = load("data/shifts_vs_eta_empirical.jld2")
#dat = d["dat"]
#N = d["N"]

## N tensor
#
# Dimensions:
# 1: Dataset index
# 2: row in the N matrix (index i), arrival state
# 3: column in the N matrix (index j), departure state
# 4: values for η = 1 /t, or ηhat, don't remember

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
    #fpath = replace(fpath, "Toninietal2016" => "Toninietal.2016")
    #fpath = replace(fpath, "Zanneetal2014rooteddated" => "Zanne.et.al.2014rooted.dated")
    println(fpath)
    tree = readtree(fpath)
    
    if all(tree.edge_length .> 0)
        bn = Base.Filesystem.basename(fpath)
        trees[bn] = tree
        datasets[bn] = SSEdata(trees[bn], ρs[bn])
    end
end

#heights = [maximum(d.node_depth) for (key, d) in datasets]

##################
#
#       FIGURE 5
#
###############

function magnitude(model, N)
    r = model.λ .- model.μ
    #shifts_weighted = N .* (model.λ .- model.λ')
    shifts_weighted = N .* (r .- r')
    mean_magnitude = sum(shifts_weighted) / sum(N)
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



n_datasets = length(datasets)
magnitudes = zeros(n_datasets)
heights = zeros(n_datasets)

for (i, name) in enumerate(keys(d))
    model = models[name]
    heights[i] = maximum(datasets[name .* ".tree"].node_depth)
    Nsum = sum(d[name]["N"], dims = 1)[1,:,:]
    m = magnitude(model, Nsum)
    magnitudes[i] = m
end

variances = zeros(n_datasets)

for (i, name) in enumerate(keys(d))
    model = models[name]
    Nsum = sum(d[name]["N"], dims = 1)[1,:,:]

    variances[i] = var_bds(model, Nsum)
end
CairoMakie.scatter(log10.(heights), log10.(variances))



fig2 = Figure(size = (350, 175), fontsize = 14, 
                figure_padding = (5,5,1,1))


xt = collect(range(extrema(magnitudes)...; length = 5))
xtl = [@sprintf("%.2f", x) for x in xt]

ax_mag_hist = Axis(fig2[1,1], 
            ylabel = L"\text{frequency}", 
            xlabel = L"\text{magnitude }(\Delta r)",
            xgridvisible = false, 
            ygridvisible = false,
            #yscale = log10, xscale = log10,
            #xticks = (xt, xtl),
            xticks = (xt, xtl),
            topspinevisible = false,
            rightspinevisible = false,
            #xticklabelrotation = π/2,
            xticklabelsize = 9,
            yticklabelsize = 9)

CairoMakie.hist!(ax_mag_hist, magnitudes, bins = 10, color = "gray")

xt = 10 .^ (collect(range(extrema(log10.(heights))...; length = 5)))
xtl = [@sprintf("%.1f", x) for x in xt]

#yt = 10 .^ (collect(range(extrema(log10.(magnitudes))...; length = 5)))
#ytl = [@sprintf("%.2f", y) for y in yt]

ax_scatter = Axis(fig2[1, 2], 
        ylabel = L"\text{magnitude } (\Delta r)", 
        #xlabel = L"\text{tree height (Ma)}",
        xgridvisible = false, 
        ygridvisible = false,
        #yscale = log10, 
        xscale = log10,
        xticks = (xt, xtl),
        #yticks = (yt, ytl),
        xlabel = L"\text{tree height (Ma)}",
        topspinevisible = false,
        rightspinevisible = false,
        #xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9)
#ylabel2 = Label(fig2[3,4], L"\text{tree height (Ma)}")



#β, Varβ, ySE = ols_regression(log10.(heights), log10.(magnitudes))
β, Varβ, ySE = ols_regression(log10.(heights), magnitudes)
linefit(x) = β[1] + β[2]*log10.(x)
x = collect(Pesto.lrange(extrema(heights)..., 20))
#yarithmetic = linefit.(log10.(x))
yarithmetic = linefit.(x)
yVar = Varβ[1,1] .+ 2 .* log10.(x) .* Varβ[1,2] .+ (log10.(x) .^ 2) .* Varβ[2,2]

#yupper = 10 .^ (yarithmetic .+ sqrt.(yVar))
#ylower = 10 .^ (yarithmetic .- sqrt.(yVar))
yupper = yarithmetic .+ 2*sqrt.(yVar)
ylower = yarithmetic .- 2*sqrt.(yVar)

#ylower[ylower .< 1.0] .= 1.0
CairoMakie.band!(ax_scatter, x, ylower, yupper, color = "#e0e0e0")

#y = 10 .^(linefit.(log10.(x)))
y = linefit.(x)
CairoMakie.lines!(ax_scatter, x, y; 
label = "OLS", markersize = 7, color = "gray", linestyle = :dash)


CairoMakie.scatter!(ax_scatter, 
                    heights, 
                    magnitudes,
                    color = "black",
                    markersize = 7)

colgap!(fig2.layout, 5)
fig2

#CairoMakie.save("figures/magnitude_empiricalbayes.pdf", fig2)
#CairoMakie.save("figures/magnitude_fixedprior.pdf", fig2)

for (i, name) in enumerate(keys(d))
    println(name, ": \t", magnitudes[i])
end

dfm = DataFrame(
    "name" => collect(keys(d)),
    "m" => magnitudes
)

sort(dfm, :m)


name = split(Base.basename(fpath), ".")[1]
x = JLD2.load("output/empirical/jld2/Sigmodontinae_VallejosGarrido2023.jld2")





