using Makie
using CSV
using DataFrames
using CairoMakie
using LaTeXStrings
using Printf

ms = ["NA", "NAN", "NULL", "NaN"]
df = CSV.read("output/rate_shift_type_munged.csv", DataFrame; missingstring = ms)
df_long = CSV.read("output/rate_shift_type_munged_long.csv", DataFrame; missingstring = ms)

## delete rows where N_sum is 0, where the fractions dont work (NA)
deleteat!(df, findall(==(0.0), df.N_sum))


df

pooled = filter(:type => x -> x == "estimate, all pooled", df)
supported = filter(:type => x -> x == "estimate, strong support", df)



#is_larger = ratios .> prior_ratios
#sum(is_larger, dims = 1)




fig = Makie.Figure(size = (450, 220), fontsize = 14, figure_padding = (1,1,1,1))

yt = collect(range(0.0, 1.0; length = 5))
ytl = [@sprintf("%.2f", y) for y in yt]

kwargs = (; 
topspinevisible = false,
rightspinevisible = false,
xgridvisible = false,
ygridvisible = false,
titlealign = :left,
xticks = (1:3, [L"\lambda", L"\mu", L"\text{both}"]),
yticks = (yt, ytl)
)


a = size(pooled)[1]
ax1 = Axis(fig[1,1]; 
    title = "a) all pooled (n = $a)",       
    ylabel = L"N_\text{rate} / N_\text{all}",
    kwargs...)

b = size(supported)[1]
ax2 = Axis(fig[1,2]; 
    title = "b) strong support (n = $b)",
    ylabel = "",
    yticklabelsvisible = false,
    kwargs...)
ylims!(ax1, (-0.1, 1.1))
ylims!(ax2, (-0.1, 1.1))

function foobar(df, ax)
    ## hard code the priors
    n = 5
    priors = [(n-1)/(n^2-1), (n-1)/(n^2-1), (n-1)^2/(n^2-1)]
    prior_ratios = zeros(size(df)[1], 3)
    for i in 1:size(prior_ratios)[1]
        prior_ratios[i,:] .= priors
    end

    arr = Matrix{Float64}(df[:,2:4])


    colors = [:steelblue, "orange", "gray"]

    xs = vcat([repeat([i], size(df)[1]) for i in 1:3]...)
    cs = vcat([repeat([i], size(df)[1]) for i in colors]...)
    
    rainclouds!(ax,
        xs .- 0.25,
        vcat(arr...); clouds = nothing,
        plot_boxplots = true,
        color = cs,
        boxplot_width = 0.5,
        side_nudge = 0.5,
        cloud_width = 0.5,
        markersize = 4,
        label = "label",
        jitter_width = 0.35)
    
    for (i, prior) in enumerate(priors)
        x = [i-0.5, i+0.5]
        y = [prior, prior]
        CairoMakie.lines!(ax, x, y, color = colors[i],
        linestyle = :solid, alpha = 0.5, label = "prior")
    end
end

foobar(pooled, ax1)
foobar(supported, ax2)

xlab = Label(fig[2,1:2], text = L"\text{type of rate shift}")

rowgap!(fig.layout, 5.0)
colgap!(fig.layout, 10.0)
fig

save("figures/rate_shift_type_simulated.pdf", fig)

