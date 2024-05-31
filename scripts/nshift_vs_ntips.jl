##################################
##
## now do it in Makie
## very annoying to do but a lot of possibility for finer control
##
##################################

function ols_regression(x, y)
    X = hcat([1 for _ in 1:length(x)], x)
    n, p = size(X)

    ## OLS
    β = (X' * X) \ X' * y
    yhat = X * β
    #sigma_squared = (1 / (n - p - 1)) * sum((y .- yhat).^2) ## MLE for sigma^2
    s_squared = (y .- yhat)' * (y .- yhat) ./ (n - p) ## OLS for sigma^2
    Varβ = inv(X' * X) .* s_squared
    yVar = x -> Varβ[1,1] + (x^2)*Varβ[2,2] + 2*x*Varβ[1,2]
    ySE = x -> sqrt.(yVar(x))
    return(β, Varβ, ySE)
end

## lrange
function lrange3(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end

using CSV
using CairoMakie

df = CSV.read("output/munged_magnitude.csv", DataFrame)


ntips = df[!,:NTips]
shifts_per_time = df[!,:N_per_time]
shifts = df[!,:N_total]


#xt = 10 .^ (collect(range(extrema(ntips)...; length = 5)))
#yt = 10 .^ (collect(range(extrema(y)...; length = 5)))
xt = collect(lrange3(extrema(Float64.(ntips))..., 5))
yt = collect(lrange3(extrema(shifts)..., 5))
yt2 = collect(lrange3(extrema(shifts_per_time)..., 5))

fmt = Printf.Format("%.0f")
xtl = [Printf.format(fmt, x) for x in xt]

fmt = Printf.Format("%.0f")
ytl = [Printf.format(fmt, y) for y in yt]

fmt = Printf.Format("%.4f")
ytl2 = [Printf.format(fmt, y) for y in yt2]

#fig1 = Figure(size=(350, 300), fontsize = 14);
fig = Figure(size = (450, 210), fontsize = 14);

## number of rate shifts (not per time)
ax1 = Axis(fig[1,1], 
            ylabel = L"\text{no. shifts }(\hat{N})", 
            xlabel = L"\text{number of tips}",
            xgridvisible = false, 
            ygridvisible = false,
            yscale = log10, xscale = log10,
            xticks = (xt, xtl),
            yticks = (yt, ytl),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelrotation = π/2,
            xticklabelsize = 9,
            yticklabelsize = 9)

CairoMakie.ylims!(ax1, minimum(yt)*0.7, maximum(yt)*1.3)

β, Varβ, ySE = ols_regression(log10.(ntips), log10.(shifts)) 
linefit(x) = β[1] + β[2]*x
x = collect(lrange3(Float64.(extrema(ntips))..., 20))
y = linefit.(log10.(x))


yupper = 10 .^ (y .+ 2 .* ySE.(log10.(x)))
ylower = 10 .^ (y .- 2 .* ySE.(log10.(x)))
ϵ = 0.7 * minimum(log10.(shifts))
ylower[ylower .< ϵ] .= ϵ

CairoMakie.band!(ax1, x, ylower, yupper, color = "#e0e0e0")
CairoMakie.scatter!(ax1, ntips, shifts; 
                    label = "asd", markersize = 7, color = "black")

x = [extrema(Float64.(ntips))...]
y = 10 .^(linefit.(log10.(x)))
CairoMakie.lines!(ax1, x, y; 
                label = "OLS", markersize = 7, color = "gray", linestyle = :dash)


## number of rate shifts per time
ax2 = Axis(fig[1,2], 
            ylabel = L"\text{shifts per time }(\hat{N}/t)", 
            xlabel = L"\text{number of tips}",
            xgridvisible = false, 
            ygridvisible = false,
            yscale = log10, xscale = log10,
            xticks = (xt, xtl),
            yticks = (yt2, ytl2),
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelrotation = π/2,
            xticklabelsize = 9,
            yticklabelsize = 9)

CairoMakie.ylims!(ax2, minimum(yt2)*0.7, maximum(yt2)*1.3)

β, Varβ, ySE = ols_regression(log10.(ntips), log10.(shifts_per_time)) 
linefit(x) = β[1] + β[2]*x
x = collect(lrange3(Float64.(extrema(ntips))..., 20))
y = linefit.(log10.(x))


yupper = 10 .^ (y .+ 2 .* ySE.(log10.(x)))
ylower = 10 .^ (y .- 2 .* ySE.(log10.(x)))
ϵ = 0.7 * minimum(log10.(shifts_per_time))
ylower[ylower .< ϵ] .= ϵ

CairoMakie.band!(ax2, x, ylower, yupper, color = "#e0e0e0")
CairoMakie.scatter!(ax2, ntips, shifts_per_time; 
                    label = "asd", markersize = 7, color = "black")

x = [extrema(Float64.(ntips))...]
y = 10 .^(linefit.(log10.(x)))
CairoMakie.lines!(ax2, x, y; 
                label = "OLS", markersize = 7, color = "gray", linestyle = :dash)


CairoMakie.save("figures/shifts_vs_ntips.pdf", fig)

fig