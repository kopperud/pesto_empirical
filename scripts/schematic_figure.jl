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


netdiv_extrema = extrema(r .- r')
model = SSEconstant(λ, μ, 1.0)

N_all = x["N"]

N_supported = N_all[is_supported,:,:]


## colors
using Colors, ColorSchemes
mycmap = ColorScheme([
    RGBf(128/255,128/255,128/255),
    RGBf(33/255,171/255,204/255),
    RGBf(204/255,84/255,0)
])

## the figure
fig = Figure(size=(650, 350), fontsize = 14, figure_padding = 0);

## the matrix
ax = Axis(fig[1,1],
        xgridvisible = false, ygridvisible = false,
        title = L"\text{b) mock rate shifts (matrix)}",
        topspinevisible = false,
        rightspinevisible = false,
        titlealign = :left,
        leftspinevisible = false,
        bottomspinevisible = false,
        xticklabelrotation = 0,
        xreversed = true
    )
hideydecorations!(ax)
hidexdecorations!(ax)
Δr_dummy = [
    NaN -0.05 -0.10 -0.15
    0.05 NaN -0.05 -0.10
    0.10 0.05 NaN -0.05
    0.15 0.10 0.05 NaN
]

heatmap!(ax, abs.(Δr_dummy); colormap=cgrad(mycmap, 3, categorical=true), aspect_ratio=1)
for i in 1:4
    for j in 1:4
        if i != j
            s = @sprintf "%1.2f" Δr_dummy[i,j]
            s = LaTeXString(s)
            text!(ax, [i-0.25],[j], text= [s], align = (:right, :center), fontsize = 18)
        end
    end
end

## the (not existing yet) mock prior
xt = [-0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15]
ax_dummy = Axis(fig[2,1],
            xgridvisible = false, ygridvisible = false,
            title = L"\text{c) mock rate shifts (frequencies)}",
            ylabel = L"\text{number of shifts}",
            xlabel = L"\text{shift size in net diversification (}\Delta r)",
            xticks = (xt, [@sprintf "%1.2f" x for x in xt]),
            titlealign = :left,
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelrotation = 0
        )

dummy_model = SSEconstant(
    [0.05, 0.10, 0.15, 0.20],
    [0.0, 0.0, 0.0, 0.0],
    0.0
)
r_dummy = dummy_model.λ
dummy_limits = extrema(r_dummy .- r_dummy')

#mids, bins = makebins(ones(4,4), dummy_model, -0.16,0.16; nbins = 6)

#offset = mids1[2] - mids[1]

#y1 = bins1[:,3] 
barplot!(ax_dummy, [-0.05, 0.05], [3.0, 3.0], color = mycmap[1], width = 0.05)
barplot!(ax_dummy, [-0.10, 0.10], [2.0, 2.0], color = mycmap[2], width = 0.05)
barplot!(ax_dummy, [-0.15, 0.15], [1.0, 1.0], color = mycmap[3], width = 0.05)

fig


ax1 = Axis(fig[1,2],
        xgridvisible = false, ygridvisible = false,
        title = L"\text{d) empirical rate shifts (all branches)}",
        ylabel = L"\text{number of shifts}",
        topspinevisible = false,
        rightspinevisible = false,
        titlealign = :left,
        xticklabelrotation = 0,
        #xlabel = L"\text{shift in net diversification (}\Delta r)",
    )
ax2 = Axis(fig[2,2],
    xgridvisible = false, ygridvisible = false,
    title = L"\text{e) empirical rate shifts (supported branches)}",
    ylabel = L"\text{number of shifts}",
    topspinevisible = false,
    titlealign = :left,
    rightspinevisible = false,
    xticklabelrotation = 0,
    xlabel = L"\text{shift size in net diversification (}\Delta r)",
)

limits = [-1.2, 1.2]

#offset = (mids[2] - mids[1])/2
w = 0.15
for (i,(N, ax)) in enumerate(zip((N_all, N_supported), (ax1, ax2)))
    NΔr = sum(N, dims = 1)[1,:,:] .* Δr
    NΔr_vec = hcat(NΔr...)
    NΔr_vec_notna = NΔr_vec[.! isnan.(NΔr_vec)]

    mids1, bins1 = makebins(ones(100,100), model, limits...; nbins = 14)
    mids2, bins2 = makebins(sum(N,dims=1)[1,:,:], model, limits...; nbins = 14)
    #offset = mids1[2] - mids[1]

    y1 = bins1[:,3] .* sum(bins2[:,3])./ sum(bins1[:,3])
    barplot!(ax, mids1, y1, color = (:gray,0.2), label = "prior", width = w)

    y2 = bins2[:,3] #./ sum(bins2[:,3])
    barplot!(ax, mids2, y2, color = (:black, 0.5), label = "posterior", width = w)

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




## the figure
fig2 = Figure(size=(350, 350), fontsize = 14)#, figure_padding = 0);


fig2











