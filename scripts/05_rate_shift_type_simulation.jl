using Revise
using BirthDeathSimulation
using ProgressMeter
using JLD2
using Distributions
using Random

function allpairwise(xs, ys)
    ny = length(xs)
    nx = length(ys)

    k = ny * nx

    λ = zeros(Base.eltype(xs), k)
    μ = zeros(Base.eltype(ys), k)

    for (i, (x, y)) in enumerate(Iterators.product(xs, ys))
        λ[i] = x
        μ[i] = y
    end

    return(λ, μ)
end


function make_quantiles(d, k)
    quantiles = zeros(k)
    step = 0.5
    for i in 1:k
        p = (i-step)/k
        quantiles[i] = Distributions.quantile(d, p)
    end
    return(quantiles)
end


#########################################
##
## Set up a 5x5 model with small rate variation
##
#########################################
λmean = 0.25
μmean = λmean * 2 / 3
sd = 0.5

λd = LogNormal(log(λmean), sd)
μd = LogNormal(log(μmean), sd)

n = 5
λq = make_quantiles(λd, n)
μq = make_quantiles(μd, n)

λ, μ = allpairwise(λq, μq)

using StatsPlots
p1 = scatter(λ, μ, xlim = (0.0, 0.4), ylim = (0.0, 0.4),
    grid = false,
    title = "model setup", label = "",
    xlabel = "speciation rate (λ)",
    ylabel = "extinction rate (μ)")

ϵ = μ ./ λ
r = λ .- μ
scatter(ϵ, r,
    grid = false,
    title = "model setup", label = "",
    xlabel = "relative extinction (μ/λ)",
    ylabel = "net-div rate (λ-μ)")

η = 0.001
model = bdsmodel(λ, μ, η)



#########################################
##
## simulate trees conditional on survival until time t
##
#########################################
Random.seed!(123)

maxtime = 100.0
maxtaxa = 10_000
n_trees = 1000
n_states = length(model.λ)
tree_heights = collect(range(30.0, 100.0; length = 8))
n_heights = length(tree_heights)

trees = Array{Tree, 1}(undef, n_trees)
starting_state = Array{Int64, 1}(undef, n_trees)
N = Array{Int64, 3}(undef, n_trees, n_states, n_states)

##################
##
## Simulation criteria (rejections)
##
##  * All lineages went extinct
##  * One of the two subtrees from the root went extinct
##  * Too many taxa (10_000)
##  * Too few taxa (100)
##  * At least one rate shift occurred
##

prog = ProgressUnknown("Complete trees sim:")
i = 1
rejection_all_extinct = 0
rejection_too_many_taxa = 0
rejection_one_subtree_extinct = 0
rejection_atleast_one_shift = 0
while i <= n_trees
    #state = rand((1:9))
    state = ((n^2)+1)÷2 ## the state with intermediate rates
    maxtime = 50.0
    tree = sim_bdshift(model, maxtime, maxtaxa, state)
    if length(tree.Leaves) < maxtaxa ## reject complete trees that termined when too many taxa
        if length(tree.Leaves) > 1 ## Reject completete trees where all taxa went extinct
            prune_extinct!(tree)
            if abs(treeheight(tree) - maxtime) < 0.001 ## reject complete trees where both root children did not survive
                N0 = +([branch.N for (idx, branch) in tree.Branches]...)

                if sum(N0) > 0
                    N[i,:,:] .= N0

                    trees[i] = tree
                    starting_state[i] = state
                    i += 1
                    ProgressMeter.next!(prog)
                else
                    rejection_atleast_one_shift += 1
                end
            else
                rejection_one_subtree_extinct += 1
            end
        else
            rejection_all_extinct += 1
        end
    else
        rejection_too_many_taxa += 1
    end
end
ProgressMeter.finish!(prog)

println("reject (too many taxa): \t", rejection_too_many_taxa)
println("reject (all taxa extinct): \t", rejection_all_extinct)
println("reject (one subtree extinct): \t", rejection_one_subtree_extinct)
println("reject (zero shifts): \t\t", rejection_atleast_one_shift)
println("accepted: \t\t\t", size(trees)[1])


y = [sum(starting_state .== i) for i in 1:9]
x = collect(range(1,9))
bar(x, y, xlabel = "starting state", ylabel = "number of trees")


heatmap(λq, μq, reshape(y, (3,3))', 
    xlabel = "speciation rate (λ)",
    ylabel = "extinction rate (μ)",
    xtick = round.(λq, digits = 3),
    ytick = round.(μq, digits = 3))
   

histogram([ntaxa(tree) for tree in trees], xlabel = "number of taxa",
        bins = 50, ylabel = "frequency (trees)")



prog = Progress(n_trees, "Writing trees:")
for i in 1:n_trees
    fpath = string(
        "data/simulations/rate_shift_type/", 
        i, ".tre"
        )
    writenewick(fpath, trees[i], model)
    next!(prog)
end
finish!(prog)


function count_shifts(N, model)
    Δλ = model.λ .- model.λ'
    Δμ = model.μ .- model.μ'

    joint_shift = (Δλ .!= 0) .& (Δμ .!= 0) 

    Δλincrease = Δλ .> 0
    Δμincrease = Δμ .> 0

    Δλdecrease = Δλ .< 0
    Δμdecrease = Δμ .< 0

    c_λ_increase = 0
    c_μ_increase = 0
    c_λ_decrease = 0
    c_μ_decrease = 0
    c_joint = 0

    K = length(model.λ)
    for i in 1:K
        for j in 1:K
            if !joint_shift[i,j]
                if Δλincrease[i,j]
                    c_λ_increase += N[i,j]
                end

                if Δμincrease[i,j]
                    c_μ_increase += N[i,j]
                end

                if Δλdecrease[i,j]
                    c_λ_decrease += N[i,j]
                end

                if Δμdecrease[i,j]
                    c_μ_decrease += N[i,j]
                end
            else
                c_joint += N[i,j]
            end
        end
    end
    return(
        [c_λ_increase,
        c_λ_decrease,
        c_μ_increase,
        c_μ_decrease,
        c_joint]
        )
end

maximum(N)
ys = [count_shifts(N[i,:,:], model) for i in 1:n_trees]

xb = [
    "λ+",
    "λ-",
    "μ+",
    "μ-",
    "λ and μ"
]

ps = []
ymax = maximum(maximum(y) for y in ys)
for i in 1:n_trees
    p = bar(xb, ys[i], 
        xlabel = "type of rate shift",  
        ylabel = "number of shifts", grid = false,
        #ylim = (0.0, ymax+1),
        title = "simulated trees",
        label = "",
        xrotation = 90)
    append!(ps, [p])
end

p5 = plot(ps[1:30]..., size = (1200, 1200))
plot(p1, ps[3])



ys2 = hcat(ys...)'

ys2
ys3 = hcat(
    sum(ys2[:,1:2], dims = 2)[:,1],
    sum(ys2[:,3:4], dims = 2)[:,1],
    ys2[:,5]
)
m = hcat([ys3[i,:] ./ sum(ys3[i,:]) for i in 1:n_trees]...)'

xls = [
    "λ", "μ", "λ and μ"
]
p = plot(grid = false, xlabel = "type of shift",
        ylabel = "frequency",
        xticks = (1:3, xls))

for i in 1:3
    boxplot!(p, [i], m[:,i], label = "", alpha = 0.5)
end
n = sqrt(n_states)
prior_single = (n-1) / (n^2-1)
plot!(p, [0.7, 1.3], [prior_single, prior_single], linewidth = 3, color = :red, linestyle = :dash,
        label = "expectation")
plot!(p, [1.7, 2.3], [prior_single, prior_single], linewidth = 3, color = :red, linestyle = :dash,
    label = "")
plot!(p, [2.7, 3.3], [1 - 2*prior_single, 1 - 2*prior_single], linewidth = 3, color = :red, linestyle = :dash,
    label = "")
p

scatter(1:300, m[:,1])
scatter!(301:600, m[:,2])
scatter!(601:900, m[:,3])

histogram(m[:,3], bins = 20)
histogram(m[:,2], bins = 20)
scatter([ntaxa(tree) for tree in trees], m[:,3])


mean(m, dims = 1)
m
r = λ .- μ
Δr = r .- r'
























































