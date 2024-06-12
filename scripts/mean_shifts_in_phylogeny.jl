using Pesto
using OrdinaryDiffEq
using Plots
using Interpolations
using Statistics

phy = readtree("/home/bkopper/projects/pesto_ms_analyses/data/primates.tre")
phy = readtree("data/empirical/Actinopterygii_Rabosky2018.tree")
#ρ = 0.635
ρ = 0.37
data = SSEdata(phy, ρ)

optres, model, i = optimize_hyperparameters(data; n = 10)
Ds, Fs = backwards_forwards_pass(model, data);


function foobar(model::SSEconstant, data::SSEdata, Ds, Fs; alg = OrdinaryDiffEq.Tsit5())
    nbranches = size(data.edges)[1]
    K = Pesto.number_of_states(model)    
    nshifts = Dict()
    ode = Pesto.shift_problem_simple(model)

    for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]
        tspan = (a,b)
        N0 = [0.0]

        p = (model.η, K, Ds[edge_idx], Fs[edge_idx]);

        prob = OrdinaryDiffEq.ODEProblem(ode, N0, tspan, p);
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = Pesto.notneg)

        nshifts[edge_idx] = sol
    end

    return(nshifts)
end


nshifts = foobar(model, data, Ds, Fs);


height = maximum(data.node_depth)+10

p1 = plot(
    xlims = (0.0, height), 
    legend = false,
    xlabel = "time before present", 
    ylabel = "N(t) per branch", 
    grid = false);
nbranches = size(data.edges)[1]
for edge_idx in 1:nbranches
    plot!(p1, nshifts[edge_idx])
end
plot!(p1, xlims = (0.0, height))
#Plots.savefig(p1, "figures/nshifts_per_branch_over_time.pdf")
Plots.savefig(p1, "figures/nshifts_per_branch_over_time.png")



ΔN_foo = Dict()

for edge_idx in 1:nbranches
    a = Ds[edge_idx].t[end]
    b = Ds[edge_idx].t[1]
    tspan = (a,b)

    Δt_target = 0.01

    number_episodes = Int64(round((a-b)/Δt_target)+1)
    Δt = (a-b) / number_episodes

    times = collect(range(a, b; length = number_episodes+1))

    ΔNs = Float64[]
    for i in 2:(number_episodes+1)
        ΔN = nshifts[edge_idx](times[i]) - nshifts[edge_idx](times[i-1])
        append!(ΔNs, ΔN)
    end

    mid_times = (times[1:end-1] .+ times[2:end]) ./ 2

    if length(mid_times) > 1
        f = linear_interpolation(reverse(mid_times), reverse(ΔNs), extrapolation_bc = Flat())
    else
        f = t -> ΔNs[1]
    end

    ΔN_foo[edge_idx] = f
end



function which_active_branches(t::Float64)
    as = [Ds[edge_idx].t[end] for edge_idx in 1:nbranches]
    bs = [Ds[edge_idx].t[1] for edge_idx in 1:nbranches]

    active_lineages = (as .> t) .& (bs .< t)

    active_branch_indices =  [idx for (idx, bool) in enumerate(active_lineages) if bool]
    return(active_branch_indices)
end


mean_ΔN_through_time = []
times2 = Float64[]

for t in range(height,0.0; length = 500)
    active_branch_indices = which_active_branches(t)
    n_active = length(active_branch_indices)

    ΔNs = Float64[]
    for edge_idx in active_branch_indices
       x = ΔN_foo[edge_idx](t) 
       append!(ΔNs, x)
    end
    
    if !isempty(ΔNs)
        mean_ΔN = Statistics.mean(ΔNs)
        append!(mean_ΔN_through_time, mean_ΔN)    
        append!(times2, t) 
    end
end


p2 = plot( grid = false,
    ylims = (0.0, 0.00004),
    ylabel = "sum(ΔN)/no. active lineages",
    xlabel = "time before the present",
    title = "fish tree (ntips = 11k)",
    xflip = true
    )
plot!(p2, times2, mean_ΔN_through_time)
Plots.savefig(p2, "figures/fish_mean_delta_N_through_time.pdf")

