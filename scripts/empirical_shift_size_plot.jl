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
#inference = "empirical"
inference = "empirical_joint"

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

meta = CSV.read("data/empirical/Phylogenies for DeepILS project.csv", DataFrame)
ρs = Dict()
for row in eachrow(meta)
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
fig = Figure(size = (480, 280), fontsize = 14, 
            figure_padding = (1,1,1,1));

colors = [:steelblue, "orange", "gray"]
labels = [
    L"\text{Shift in \lambda}", 
    L"\text{Shift in \mu}", 
    L"\text{Shift in both}"
    ]

g = fig[1,1] = GridLayout()

name_subset = [
    "Actinopterygii_Rabosky2018",
    "Chondrichthyes_Stein2018",
    "Primates_Springer2012",
    "Sigmodontinae_VallejosGarrido2023",
    "Mammalia_AlvarezCarretero2022",
    "Rosidae_Sun2020",
    #"Rhopalocera_Kawahara2023",
    #"Monocots_Howard2019",
    "Polypodiophyta_Nitta2022",
    "Agaricomycetes_SanchezGarcia2020",
]

titles = []
dnames = collect(keys(datasets))
for dname in dnames
    append!(titles, [split(dname, "_")[1]])
end


#datasets_order = [34, 9, 28, 29, 4]

axs = []
q = 1
for i in 1:3, j in 1:3
    if (i == 1) & (j == 3)
        #ax = Axis(fig[i,j], scene = false)
        #hidedecorations!(ax)
    else
        if i < 3
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

        if i < 3
            hidexdecorations!(ax, ticks = false)
        end
    
        #=if j > 1
            hideydecorations!(ax, ticks = false, ticklabels = false)
        end=#

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
    nbins = 20

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
leg = Legend(g[1,3], elements, labels, labelsize=9, 
            patchsize = (10.0f0, 10.0f0),
            framevisible = false)

leg.tellheight = true
leg.tellwidth = true


linkxaxes!(axs...)
ylabel = Label(g[1:3, 0], L"\text{number of rate shifts }(\hat{N})", rotation = π/2)
xlabel = Label(g[4, 1:3], L"\text{shift size in net diversification }(\Delta r)")




for col in 1:3
    colsize!(g, col, Relative(0.30))
end
colsize!(g, 0, Relative(0.1))

for i in 1:3
    rowsize!(g, i, Relative(0.30))
end
rowsize!(g, 4, Relative(0.1))

colgap!(g, 5)
rowgap!(g, 3)


#linkaxes!(ax_scatter1, ax_scatter2)

fig
#set_theme!(fig, figure_padding = 0)
#CairoMakie.save("figures/fig1_fixedprior.pdf", fig)
#CairoMakie.save("figures/fig1_empiricalbayes.pdf", fig)
#CairoMakie.save("figures/fig1_empirical_joint.pdf", fig)
#CairoMakie.save("figures/histogram_magnitude_etaoneshift.pdf", fig)




###
## plot each phylogeny individually
##



for name in keys(models)  
    model = models[name]

    #Nmatrix = N[dataset_index,:,:,1]
    Nmatrix = sum(d[name]["N"], dims = 1)[1,:,:]
    nbins = 20

    Ns = zeros(3,nbins)
    filters = ["extinction", "speciation", ""]

    dfs = []
    #netdiv_extrema = [-1.2, 1.2]
    r = model.λ .- model.μ
    Δr = r .- r'
    netdiv_extrema = extrema(Δr)

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
    #ratios[q,:] .= rs

    bar_df = vcat(dfs...)

    fig = Figure(size = (200,200))

    ax = Axis(fig[1,1], 
            xgridvisible = false, 
            ygridvisible = false,
            title = name,
            titlesize = 9,
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelrotation = π/2,
            xticklabelsize = 9,
            xlabel = "shift size (netdiv)",
            ylabel = "number of shifts",
            yticklabelsize = 9)

    barplot!(ax, 
        bar_df[!, "mids"], 
        bar_df[!, "Δr"], 
        stack = bar_df[!, "subset"],
        color = [colors[x] for x in bar_df[!, "subset"]],
        )

    
    save("/tmp/empirical_shiftfigs/$name.pdf", fig)    

end


###
## plot them in a big plot
##

fig3 = Figure(size = (800, 800));

q = 1
qnames = [keys(models)...]
axs = []
for i in 1:7
    for j in 1:7
        if q <= length(keys(models))
            ax = Axis(fig3[i,j], 
                xgridvisible = false, 
                ygridvisible = false,
                titlesize = 7,
                title = qnames[q],
                topspinevisible = false,
                rightspinevisible = false,
                xticklabelrotation = π/2,
                xticklabelsize = 9,
                yticklabelsize = 9)
            push!(axs, ax)
            q += 1
        end
    end
end

i = 1
for name in keys(models)
    model = models[name]
    Nmatrix = sum(d[name]["N"], dims = 1)[1,:,:]
    nbins = 20

    Ns = zeros(3,nbins)
    filters = ["extinction", "speciation", ""]

    dfs = []
    if true
        netdiv_extrema = [-1.5, 1.5]
    else
        r = model.λ .- model.μ
        Δr = r .- r'
        netdiv_extrema = extrema(Δr)
    end

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
    #ratios[q,:] .= rs

    bar_df = vcat(dfs...)


    ax = axs[i]
    barplot!(ax, 
        bar_df[!, "mids"], 
        bar_df[!, "Δr"], 
        stack = bar_df[!, "subset"],
        color = [colors[x] for x in bar_df[!, "subset"]],
        )
       
    i += 1
end

#CairoMakie.save("figures/shiftsize_free_xlimit.pdf", fig3)
CairoMakie.save("figures/shiftsize_fixed_xlimit.pdf", fig3)




## compute the variance
kingdom = Dict{String, String}()
for row in eachrow(meta)
    fn = row[:Filename]
    if ismissing(fn)
        continue
    end
    name = replace(fn, ".tree" => "")
    kingdom[name] = row[:Kingdom]
end



name = "Sigmodontinae_VallejosGarrido2023"


mags = Float64[]
vars = Float64[]
heights = Float64[]
names = String[]
kingdoms = String[]
netdivs = Float64[]
sp_rates = Float64[]

@showprogress for name in keys(models)

    model = models[name]
    N = sum(d[name]["N"], dims = 1)[1,:,:]

    Δr = model.λ .- model.μ

    mag = sum(Δr .* N) / sum(N)
    var = sum((Δr .* N .- mag) .^2) #/ sum(N)
    push!(mags, mag)
    push!(vars, var)
    height = maximum(datasets[name * ".tree"].node_depth)
    push!(heights, height)
    push!(names, name)
    push!(kingdoms, kingdom[name])

    fname = "output/" * inference * "/newick/" * name * ".tre"
    @rput fname
    R"""
    library(treeio)
    tr <- read.beast.newick(fname)
    xdf <- tr@data
    xdf <- xdf[order(xdf$edge),]
    netdiv <- sum(xdf$mean_netdiv * tr@phylo$edge.length) / sum(tr@phylo$edge.length)
    speciationrate <- sum(xdf$mean_lambda * tr@phylo$edge.length) / sum(tr@phylo$edge.length)
    """
    @rget netdiv
    @rget speciationrate

    push!(netdivs, netdiv)
    push!(sp_rates, speciationrate)
end

newdf = DataFrame(
    "magnitude" => mags,
    "variance" => vars,
    "height" => heights,
    "name" => names,
    "kingdom" => kingdoms,
    "mean_netdiv" => netdivs,
    "mean_speciation" => sp_rates,
)

plants = filter(:kingdom => x -> x == "Plantae", newdf)
animals = filter(:kingdom => x -> x == "Animalia", newdf)
fungi = filter(:kingdom => x -> x == "Fungi", newdf)

fig2 = Figure();
ax1 = Axis(fig2[1,1], 
        xscale = log10, 
        yscale = log10,
        xgridvisible = false, 
        ygridvisible = false,
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xlabel = L"\text{tree height (Ma)}",
        ylabel = L"\text{variation in shift size}",
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9);
ax2 = Axis(fig2[1,2], 
        xscale = log10, 
        yscale = log10,
        xgridvisible = false, 
        ygridvisible = false,
        titlesize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        xlabel = L"\text{mean sp rate}",
        ylabel = L"\text{variation in shift size}",
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9);
#scatter!(ax, heights, sqrt.(vars))
for (kng, label) in zip((plants, animals, fungi), ("Plantae", "Animalia", "Fungi"))
    scatter!(ax1, kng[!,:height], kng[!,:variance], label = label)
    scatter!(ax2, kng[!,:mean_speciation], kng[!,:variance], label = label)
end
axislegend(ax1; position = :lt);
#axislegend(ax2; position = :ct);
fig2



