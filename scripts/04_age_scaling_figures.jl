using Makie
using CSV
using DataFrames
using CairoMakie
using LaTeXStrings
using Printf

ms = ["NA", "NAN", "NULL", "NaN"]
df = CSV.read("output/age_scaling_effect_munged.csv", DataFrame; missingstring = ms)

#df = df[df[!,:type] .== "pooled",:]
df = df[df[!,:type] .== "pooled",:]
#df = df[df[!,:inference] .== "age_scaling_effect",:]
df = df[df[!,:inference] .== "age_scaling_effect",:]

## filter for at least one strongly supported shift
#df = df[df[!,:how_many_supported ] .> 0,:]

df[:,["N_total", "N_per_time", "how_many_supported", "support_per_time"]]

df_true = CSV.read("output/age_scaling_effect_true_nshift.csv", DataFrame)
df_true[!,:height] = Int64.(round.(df_true[!,:height]))


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


#### TODO
## need to replicate the scatter plot regression, 
## and calculate the slope for the simulated data


fig3 = Figure(size = (650, 300))

heights = [30, 40, 50, 60, 70, 80, 90, 100]
xt = heights
ax = Axis(
    fig3[1,1],
    xticks = xt,
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    yscale = log10,
    xlabel = L"\text{tree height (Ma)}",
    ylabel = L"\text{no. rate shift events per time }(N/t)",
)

ylims!(ax, 1e-6, 0.2)


colors = [:black, :orange]
labels = ["estimated shifts", "true shifts"]

for (i, this_df) in enumerate([df, df_true])
    for height in heights
        df1 = filter(:height => h -> h == height, this_df)

        y = df1[!,:N_per_time]
        x = df1[!,:height]
        w = 3.0
        x_jitter = x .+ (rand(500).-0.5) .* w
        #x_jitter 
        if i == 1
            x_jitter .-= w/2.0
        else
            x_jitter .+= w/2.0
        end

        scatter!(ax, x_jitter, y, color = colors[i], label = labels[i], markersize = 4)
    end
end

η_true = 0.0008
lines!(ax, [25, 105], [η_true, η_true], color = :red, linestyle = :dash)

elem_1 = MarkerElement(color = :orange, marker = :circle, markersize = 15)
elem_2 = MarkerElement(color = :black, marker = :circle, markersize = 15)
elem_3 = LineElement(color = :red, linestyle = :dash)

Legend(fig3[1,2], 
    [elem_1, elem_2, elem_3],
    ["true shifts", "estimated shifts", "shift rate (true)"],
    patchsize = (35, 35), rowgap = 5,
    framevisible = false
    )
   
fig3
save("figures/age_scaling_effect2.pdf", fig3)













































Makie.hist(df[df[!,:height] .== 30,:N_per_time])

Makie.scatter(df[!,:height], df[!,:N_per_time])










