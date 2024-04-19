using Makie
using CSV
using DataFrames
using CairoMakie
using LaTeXStrings
using Printf
using Statistics

ms = ["NA", "NAN", "NULL", "NaN"]
df = CSV.read("output/rate_shift_type_munged.csv", DataFrame; missingstring = ms)
df_long = CSV.read("output/rate_shift_type_munged_long.csv", DataFrame; missingstring = ms)

## delete rows where N_sum is 0, where the fractions dont work (NA)
deleteat!(df, findall(==(0.0), df.N_sum))

df

pooled = filter(:type => x -> x == "estimate, all pooled", df)
#supported = filter(:type => x -> x == "estimate, strong support", df)




fig = Makie.Figure(size = (400, 250), fontsize = 14)
#figure_padding = (1,1,1,1))

yt = collect(range(0.0, 1.0; length = 5))
ytl = [@sprintf("%.2f", y) for y in yt]

kwargs = (; 
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    #titlealign = :left,
    xticks = (1:3, [L"\lambda", L"\mu", L"\text{both}"]),
    yticks = (yt, ytl)
)


a = size(pooled)[1]
ax1 = Axis(fig[1,1]; 
    #title = "a) all pooled (n = $a)",       
    ylabel = L"N_\text{rate} / N_\text{all}",
    xlabel = L"\text{type of rate shift}",
    kwargs...)

#=b = size(supported)[1]
ax2 = Axis(fig[1,2]; 
    title = "b) strong support (n = $b)",
    ylabel = "",
    yticklabelsvisible = false,
    kwargs...) =#
ylims!(ax1, (-0.1, 1.1))
#ylims!(ax2, (-0.1, 1.1))

function foobar(df, ax)
    ## hard code the priors
    n = 6
    priors = [(n-1)/(n^2-1), (n-1)/(n^2-1), (n-1)^2/(n^2-1)]
    prior_ratios = zeros(size(df)[1], 3)
    for i in 1:size(prior_ratios)[1]
        prior_ratios[i,:] .= priors
    end

    arr = Matrix{Float64}(df[:,["N_lambda_ratio", "N_mu_ratio", "N_both_ratio"]])


    colors = [:steelblue, "orange", "gray"]

    xs = vcat([repeat([i], size(arr)[1]) for i in 1:3]...)
    cs = vcat([repeat([i], size(arr)[1]) for i in colors]...)
    
    rainclouds!(ax,
        xs .- 0.3,
        vcat(arr...); clouds = nothing,
        plot_boxplots = true,
        color = cs,
        boxplot_width = 0.5,
        side_nudge = 0.6,
        cloud_width = 0.6,
        markersize = 4,
        label = "label",
        jitter_width = 0.40)
    
    for (i, prior) in enumerate(priors)
        x = [i-0.5, i+0.5]
        y = [prior, prior]
        CairoMakie.lines!(ax, x, y, color = colors[i],
        linestyle = :solid, alpha = 0.5, label = "prior")
    end
end

foobar(pooled, ax1)
#foobar(supported, ax2)

#xlab = Label(fig[2,1], text = L"\text{type of rate shift}")

#rowgap!(fig.layout, 5.0)
#colgap!(fig.layout, 10.0)
fig

#save("figures/rate_shift_type_simulated.pdf", fig)
save("figures/rate_shift_type_1.pdf", fig)


hist(df[!,"N_mu"])

#(df, "output/rate_shift_type_true_shifts.csv")



df_true = CSV.read("output/rate_shift_type_true_shifts.csv", DataFrame; missingstring = ms)

df_true[!,"N_lambda_true"] = df_true[!,"N_lambda"]
df_true[!,"N_mu_true"] = df_true[!,"N_mu"]
df_true[!,"N_both_true"] = df_true[!,"N_both"]

df_true[!,"N_lambda_ratio_true"] = df_true[!,"N_lambda_ratio"]
df_true[!,"N_mu_ratio_true"] = df_true[!,"N_mu_ratio"]
df_true[!,"N_both_ratio_true"] = df_true[!,"N_both_ratio"]


df_true[!,"N_sum_true"] = df_true[!,"N_sum"]

#hcat(df_true, df; makeunique=true)
dfx = innerjoin(df, df_true, on = :tree_index; makeunique = true)



function lrange2(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end

mean(dfx[!,:N_lambda] .- dfx[!,:N_lambda_true])
ntip = dfx[!,:ntip]

xt = collect(lrange2(10.0, 10_000.0, 4))
xtl = [@sprintf("%.0f", x) for x in xt]

kwargs2 = (; 
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
)


error = dfx[!,:N_lambda] .- dfx[!,:N_lambda_true]
relative_error = error ./ dfx[!,:N_sum_true]


#################
fig2 = Figure(size = (650, 400))

ax1 = Axis(fig2[1,1];
        title = L"\text{true}", 
        kwargs2...)
ax2 = Axis(fig2[2,1];
        title = L"\text{estimated}",
        xlabel = L"\text{number of trees}",
        kwargs2...)
ax3 = Axis(fig2[1,3];
        title = L"\text{true}",
        kwargs2...)
ax4 = Axis(fig2[2,3];
        title = L"\text{estimated}",
        xlabel = L"\text{number of trees}",
        kwargs2...)

N_rate = dfx[!,:N_lambda]
N_rate_true = dfx[!,:N_lambda_true]
h_sp = hist!(ax1, N_rate_true, label = "true", bins = 20, direction = :x)
hist!(ax2, N_rate, label = "estimated", bins = 30, direction = :x)

N_rate = dfx[!,:N_mu]
N_rate_true = dfx[!,:N_mu_true]
h_mu = hist!(ax3, N_rate_true, label = "true", bins = 20, color = :orange, direction = :x)
hist!(ax4, N_rate, label = "estimated", bins = 30, color = :orange, direction = :x)

linkxaxes!(ax1, ax2, ax3, ax4)
linkyaxes!(ax1, ax3)
linkyaxes!(ax2, ax4)

Label(fig2[1:2,0], L"\text{no. shifts in speciation rate (}N_\lambda)", rotation = pi/2)
Label(fig2[1:2,2], L"\text{no. shifts in extinction rate (}N_\mu)", rotation = pi/2)

Legend(fig2[1:2,4], 
    [h_sp, h_mu], 
    [L"\text{speciation}", L"\text{extinction}"]
    )
fig2
save("figures/rate_shift_type_2.pdf", fig2)


### now for the ratios
error_lambda = dfx[!,:N_lambda_ratio] .- dfx[!,:N_lambda_ratio_true]
error_mu = dfx[!,:N_mu_ratio] .- dfx[!,:N_mu_ratio_true]


fig3 = Figure(size = (600, 400))
ax1 = Axis(fig3[1,1];
    xlabel = L"\text{number of trees}",
    title = L"\text{shifts in speciation rate}",
    kwargs2...)
ax2 = Axis(fig3[1,2];
    xlabel = L"\text{number of trees}",
    title = L"\text{shifts in extinction rate}",
    yticklabelsvisible = false,
    kwargs2...)
ax3 = Axis(fig3[2,1];
    xlabel = L"\text{number of tips}",
    xt = [0, ]
    kwargs2...)
ax4 = Axis(fig3[2,2];
    xlabel = L"\text{number of tips}",
    yticklabelsvisible = false,
    kwargs2...)

h_sp = hist!(ax1, error_lambda, bins = 30, color = :steelblue, label = "estimation error", direction = :x)
line = lines!(ax1, [0.0, 400], [0.0, 0.0], linestyle = :dash, color = :black, label = "zero error")

h_mu = hist!(ax2, error_mu, bins = 30, color = :orange, label = "estimation error", direction = :x)
lines!(ax2, [0.0, 400], [0.0, 0.0], linestyle = :dash, color = :black, label = "zero error")
prior_ratio = (4/24) .- (5/35)

linkaxes!(ax1, ax2)
linkaxes!(ax3, ax4)

scatter!(ax3, ntip, error_lambda, markersize = 5, color = :steelblue, label = "estimation error")
lines!(ax3, [extrema(ntip)...], [0.0, 0.0], linestyle = :dash, color = :black, label = "zero error")

scatter!(ax4, ntip, error_mu, markersize = 5, color = :orange, label = "estimation error")
lines!(ax4, [extrema(ntip)...], [0.0, 0.0], linestyle = :dash, color = :black, label = "zero error")
Label(fig3[1:2,0],  L"\text{estimation error }\left (\frac{\hat{N}_\lambda}{\hat{N}} - \frac{N_{\lambda,\text{true}}}{N_\text{true}} \right )", rotation = pi/2)
Legend(fig3[1:2,3], 
    [h_sp, h_mu, line], 
    [L"\text{speciation}", L"\text{extinction}", L"\text{zero error}"]
    )

rowgap!(fig3.layout, 5.0)
colgap!(fig3.layout, 10.0)
fig3

save("figures/rate_shift_type_3.pdf", fig3)



fig5 = Figure()
ax = Axis(fig5[1,1])
hist!(ax, dfx[!,:N_lambda_ratio_true], bins = 30)
hist!(ax, dfx[!,:N_lambda_ratio], bins = 8)
fig5







