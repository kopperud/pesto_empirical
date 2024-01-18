using Makie
using CSV
using DataFrames
using CairoMakie
using LaTeXStrings
using Printf

ms = ["NA", "NAN", "NULL", "NaN"]
df = CSV.read("output/age_scaling_effect_munged.csv", DataFrame; missingstring = ms)

df

df[:,["N_total", "N_per_time", "how_many_supported", "support_per_time"]]


fig = Makie.Figure(size = (450, 220), fontsize = 14, figure_padding = (1,1,1,1))

xt = range(30, 100; length = 8) |> collect

kwargs = (; 
topspinevisible = false,
rightspinevisible = false,
xgridvisible = false,
ygridvisible = false,
titlealign = :left,
#yscale = Makie.sqrt,
#xticks = (1:3, [L"\lambda", L"\mu", L"\text{both}"]),
xticks = (xt, [@sprintf("%.0f", x) for x in xt])
#yticks = (yt, ytl)
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


fig2 = Makie.Figure(size = (450, 220), fontsize = 14, figure_padding = (1,1,1,1))

ax2 = Axis(fig2[1,1]; 
    ylabel = L"N/t",
    xlabel = "tree height (Ma)",
    #title = "a) all pooled (n = $a)",       
    #ylabel = L"N_\text{rate} / N_\text{all}",
    kwargs...)

for h in unique(df[!,:height])
    y = df[df[!,:height] .== h,:N_per_time]
    #hist!(ax2, y, offset = h, scale_to = -5, direction = :x)
    boxplot!(repeat([h], 1000), y, show_outliers = false, width = 5)
    #boxplot!(repeat([h+2.5], 1000), y, show_outliers = true, width = 5)
end

fig2


function foobaz(x::T, y) where {T<:Real}
    x+y
end

function foobaz(x::T, y::T, z) where {T<:Real}
    x+y+z
end

@code_warntype foobaz(1.0, 2.0, 3.0)

@code_lowered foobaz(1.0, 2.0)

@code_native foobaz(1.0, 2.0)
@code_native foobaz(Float16(1.0), Float16(2.0))

l = Float64[]
append!(l, Int64(0))


function barfoo()
    r = rand(1)[1]

    if r > 0.5
        x = 3.0
    else
        x = [3.0, 2.0]
    end

    return(x)
end

@code_warntype barfoo()






































Makie.hist(df[df[!,:height] .== 30,:N_per_time])

Makie.scatter(df[!,:height], df[!,:N_per_time])










