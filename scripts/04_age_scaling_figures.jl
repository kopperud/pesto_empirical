using Makie
using CSV
using DataFrames
using CairoMakie
using LaTeXStrings
using Printf

ms = ["NA", "NAN", "NULL", "NaN"]
df = CSV.read("output/age_scaling_effect_munged.csv", DataFrame; missingstring = ms)

#df = df[df[!,:type] .== "pooled",:]
df = df[df[!,:type] .== "strong support",:]

## filter for at least one strongly supported shift
df = df[df[!,:how_many_supported ] .> 0,:]

df[:,["N_total", "N_per_time", "how_many_supported", "support_per_time"]]


fig = Makie.Figure(size = (450, 220), fontsize = 14, figure_padding = (1,1,1,1))

xt = range(30, 100; length = 8) |> collect

kwargs = (; 
topspinevisible = false,
rightspinevisible = false,
xgridvisible = false,
ygridvisible = false,
titlealign = :center,
xticks = (xt, [@sprintf("%.0f", x) for x in xt])
)

ax1 = Axis(fig[1,1]; 
    ylabel = L"N/t",
    xlabel = "tree height (Ma)",
    #title = "a) all pooled (n = $a)",       
    #ylabel = L"N_\text{rate} / N_\text{all}",
    kwargs...)

#xs = 30

colors = [:steelblue, "orange", "gray"]
#cs = vcat([repeat([i], size(df)[1]) for i in colors]...)

rainclouds!(ax1, 
    df[!,:height] .- 2,
    df[!,:N_per_time]; clouds = nothing,
    plot_boxplots = true,
    #color = cs,
    boxplot_width = 4,
    side_nudge = 4,
    cloud_width = 4,
    markersize = 4,
    label = "label",
    jitter_width = 4)

fig


fig2 = Makie.Figure(size = (450, 450), fontsize = 14, figure_padding = (1,1,1,1))

n = size(df)[1]
ax1 = Axis(fig2[1,1]; 
    ylabel = L"\hat{N}/t \text{ (all branches)}",
    #xlabel = "tree height (Ma)",
    title = L"\text{simulated trees, filtered for }\geq 1\text{ supported branch }(n = %$n)",
    titlealign = :center,
    kwargs...)

ax2 = Axis(fig2[2,1]; 
    ylabel = L"\text{no. supported branches}/t",
    xlabel = L"\text{tree height (Ma)}",
    kwargs...)

for h in unique(df[!,:height])
    y = df[df[!,:height] .== h,:N_per_time]
    hist!(ax1, y, offset = h, scale_to = -5, direction = :x, bins = 50)
    boxplot!(ax1, repeat([h+2.5], length(y)), y, show_outliers = true, width = 5)

    y = df[df[!,:height] .== h,:support_per_time]
    hist!(ax2, y, offset = h, scale_to = -5, direction = :x, bins = 50)
    boxplot!(ax2, repeat([h+2.5], length(y)), y, show_outliers = true, width = 5)
    #boxplot!(repeat([h+2.5], 1000), y, show_outliers = true, width = 5)
end

#ylims!(ax2, -0.001, 0.03)
rowgap!(fig2.layout, 1.0)
fig2

save("figures/age_scaling_effect.pdf", fig2)












































Makie.hist(df[df[!,:height] .== 30,:N_per_time])

Makie.scatter(df[!,:height], df[!,:N_per_time])










