using Distributions
using CairoMakie
using LaTeXStrings
using Pesto

## the priors
λml = 0.15
μml = 0.10

H = 0.587
n = 10

dλ = LogNormal(log(λml), 2*H)
dμ = LogNormal(log(µml), 2*H)

λquantiles = make_quantiles(dλ, n)
µquantiles = make_quantiles(dμ, n)
λ, μ = allpairwise(λquantiles, µquantiles)

r = λ .- μ

Δr = r .- r'

for i in 1:(n^2)
    Δr[i,i] = NaN
end


x = load("output/empirical_fixedprior/jld2/Actinopterygii_Rabosky2018.jld2")
df = CSV.read("output/empirical_fixedprior/rates/Actinopterygii_Rabosky2018.csv", DataFrame)
sorted_df = sort(df, :edge)
is_supported1 = (sorted_df[!,:shift_bf] .> 10) .& (sorted_df[!,:nshift] .> 0.5)
is_supported = is_supported1[2:end]


netdiv_extrema = extrema(r)
model = SSEconstant(λ, μ, 1.0)

N_all = x["N"]

N_supported = N_all[is_supported,:,:]




## the figure
fig = Figure(size=(350, 350), fontsize = 14, figure_padding = 0);
ax1 = Axis(fig[1,1],
        xgridvisible = false, ygridvisible = false,
        title = L"\text{d) empirical prior (all branches)}",
        ylabel = L"\text{number of shifts}",
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelrotation = 0,
        #xlabel = L"\text{shift in net diversification (}\Delta r)",
    )
ax2 = Axis(fig[2,1],
    xgridvisible = false, ygridvisible = false,
    title = L"\text{e) empirical prior (supported branches)}",
    ylabel = L"\text{number of shifts}",
    topspinevisible = false,
    rightspinevisible = false,
    xticklabelrotation = 0,
    xlabel = L"\text{shift in net diversification (}\Delta r)",
)


for (i,(N, ax)) in enumerate(zip((N_all, N_supported), (ax1, ax2)))
    NΔr = sum(N, dims = 1)[1,:,:] .* Δr
    NΔr_vec = hcat(NΔr...)
    NΔr_vec_notna = NΔr_vec[.! isnan.(NΔr_vec)]

    mids1, bins1 = makebins(ones(100,100), model, netdiv_extrema...; nbins = 14)
    mids2, bins2 = makebins(sum(N,dims=1)[1,:,:], model, netdiv_extrema...; nbins = 14)
    #offset = mids1[2] - mids[1]

    y1 = bins1[:,3] .* sum(bins2[:,3])./ sum(bins1[:,3])
    barplot!(ax, mids1, y1, color = (:gray,0.2), label = "prior")

    y2 = bins2[:,3] #./ sum(bins2[:,3])
    barplot!(ax, mids2, y2, color = (:red, 0.5), label = "posterior")

    Δr_vec = hcat(Δr...)
    Δr_vec_notna = Δr_vec[.! isnan.(Δr_vec)]

end
#xlabel = Label(fig1[0,1:2], L"\text{shift in net diversification (}\Delta r)")

#fig[1,1] = Legend(fig, ax, "", framevisible = false, tellwidth = false, tellheight = true, position = :rt)
axislegend(ax1, framevisible = false)
fig


CairoMakie.save("figures/empirical_shiftprior.pdf", fig)

using JLD2

extrema(NΔr_vec_notna)

hist(NΔr_vec_notna)


bins





sum(N)
sum(NΔr_vec_notna)



fig













