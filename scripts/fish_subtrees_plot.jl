using Glob
using DataFrames
using CSV
using JLD2
using CairoMakie
using LaTeXStrings
using Printf

fpaths = Glob.glob("output/empirical_fish_subtrees/rates/*.csv")

dfs = DataFrame[]

for fpath in fpaths
    df = CSV.read(fpath, DataFrame)

    node = parse(Int64, split(replace(Base.basename(fpath), ".csv" => ""), "_")[end])

    df[!,:node] = repeat([node], size(df)[1])
    push!(dfs, df)
end

fpaths = Glob.glob("output/empirical_fish_subtrees/jld2/*.jld2")


treeheights = [load(fpath, "treeheight") for fpath in fpaths]
treelengths = [load(fpath, "treelength") for fpath in fpaths]
ntips = [load(fpath, "ntip") for fpath in fpaths]
N = [sum(df[!,:nshift]) for df in dfs]

N_per_time = N ./ treelengths



function lrange(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end


xt = lrange(2.0, 300.0, 8)
xtl = [@sprintf "%2.1f" x for x in xt]
ms = 5

fig = Figure(size = (500, 500))
ax1 = Axis(fig[1,1],
    xscale = log10,
    yscale = log10,
    xticks = (xt, xtl),
    xticklabelrotation = π/2,
    ygridvisible = false,
    xgridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
    title = L"\text{a) full y}\endash \text{axis}",
    titlealign = :left,
    ylabel = L"\hat{N}/t",
    xlabel = L"\text{tree height (Ma)}",
)

ax2 = Axis(fig[1,2],
    xscale = log10,
    yscale = log10,
    xticks = (xt, xtl),
    xticklabelrotation = π/2,
    ygridvisible = false,
    xgridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
    title = L"\text{b) cut off y}\endash \text{axis}",
    titlealign = :left,
    #title = "cut off y-axis",
    #ylabel = L"\hat{N}/t",
    xlabel = L"\text{tree height (Ma)}",
)
ylims!(ax2, (0.7*1e-3, 2.0))

scatter!(ax1, treeheights, N_per_time, color = :black, markersize = ms)
scatter!(ax2, treeheights, N_per_time, color = :black, markersize = ms)

xt = lrange(30.0, 10000.0, 8)
xtl = [@sprintf "%2.1f" x for x in xt]


ax3 = Axis(fig[2,1],
    xscale = log10,
    yscale = log10,
    xticks = (xt, xtl),
    xticklabelrotation = π/2,
    ygridvisible = false,
    xgridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
    title = L"\text{c)}",
    titlealign = :left,
    #title = "full y-axis",
    ylabel = L"\hat{N}/t",
    xlabel = L"\text{number of tips}",
)
scatter!(ax3, ntips, N_per_time, color = :black, markersize = ms)


ax4 = Axis(fig[2,2],
    xscale = log10,
    yscale = log10,
    xticks = (xt, xtl),
    xticklabelrotation = π/2,
    ygridvisible = false,
    xgridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
    title = L"\text{d)}",
    titlealign = :left,
    #title = "cut off y-axis",
    #ylabel = L"\hat{N}/t",
    xlabel = L"\text{number of tips}",
)
ylims!(0.7*1e-3, 2.0)
scatter!(ax4, ntips, N_per_time, color = :black, markersize = ms)

rowgap!(fig.layout, 7)
colgap!(fig.layout, 7)

fig

save("figures/fishtree_subclades.pdf", fig)




fig2 = Figure()
yt = lrange(1.0, 300.0, 6)
ytl = [@sprintf "%2.1f" y for y in yt]
xt = lrange(1.0, 300.0, 6)
xtl = [@sprintf "%2.1f" x for x in xt]

ax = Axis(fig2[1,1], 
    xscale = log10, 
    yscale = log10,
    ygridvisible = false,
    xgridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
    yticks = (yt, ytl),
    xticks = (xt, xtl),
    ylabel = L"\text{\hat{N}}",
    xlabel = L"\text{tree height (Ma)}",
    )
CairoMakie.scatter!(ax, treeheights, N,
    color = :black
)
ylims!(ax, 0.5, 500)
fig2





CairoMakie.hist(N[N .< 1.0])




argmin((N_per_time .- 1e-03) .^ 2)






