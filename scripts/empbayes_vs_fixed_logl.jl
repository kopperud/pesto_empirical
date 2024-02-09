

using Glob
using JLD2

using CairoMakie

fpaths = Glob


logls = Dict()
etas = Dict()

for inference in ["empirical", "empirical_fixedprior"]
    fpaths = Glob.glob(inference * "/jld2/*.jld2", "/home/bkopper/hpc_pesto_empirical/")
  
    @showprogress for fpath in fpaths
        name = split(Base.basename(fpath), ".")[1]
        logl = JLD2.load(fpath, "logl")
        logls[name, inference] = logl
        eta = JLD2.load(fpath, "etaml")
        etas[name, inference] = eta
    end
end



names = [split(Base.basename(fpath), ".")[1] for fpath in fpaths]

x1 = [logls[name, "empirical"] for name in names]
x2 = [logls[name, "empirical_fixedprior"] for name in names]
eta1 = [etas[name, "empirical"] for name in names]
eta2 = [etas[name, "empirical_fixedprior"] for name in names]
Δlogl = x1 .- x2
x = [i for i in 1:length(x2)]

heights = [maximum(datasets[name * ".tree"].node_depth) for name in names]
sampling_fractions = [maximum(datasets[name * ".tree"].ρ) for name in names]

fig = CairoMakie.Figure(size=(800,400))
ax = Axis(fig[1,1], 
    ylabel = "logl (emp Bayes) - logl (fixed)",
)
ylims!(ax, -500, 500)
scatter!(ax, sampling_fractions, Δlogl)

ax2 = Axis(fig[1,2], 
    ylabel = "logl (emp Bayes) - logl (fixed)",
    xscale = log10
)
scatter!(ax2, eta1, Δlogl)

ax3 = Axis(fig[1,3], 
    ylabel = "logl (emp Bayes) - logl (fixed)",
    xscale = log10
)
scatter!(ax3, heights, eta1 .- eta2)

fig

names[argmin(Δlogl)]

#JLD2.load(fpaths[end-1], "logl")

fpaths[end]






