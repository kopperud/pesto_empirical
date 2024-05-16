using Revise
using BirthDeathSimulation
using Distributions
using Random

## target metrics:
# ntips : 11637
# N = 495
# age: 368
# lambda: 0.0927
# mu: 0.0102

function make_quantiles(d, k)
    quantiles = zeros(k)
    step = 0.5
    for i in 1:k
        p = (i-step)/k
        quantiles[i] = Distributions.quantile(d, p)
    end
    return(quantiles)
end

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



Random.seed!(12345)

λ = [0.05, 0.08, 0.06]
μ = [0.03, 0.06, 0.03]

sd = 0.6
dλ = LogNormal(log(0.0927), sd)
dμ = LogNormal(log(0.0102), sd)

n = 5
λquantiles = make_quantiles(dλ, n)
μquantiles = make_quantiles(dμ, n)

λ,μ = allpairwise(λquantiles, μquantiles)





η = 0.002179824
bds3 = bdsmodel(λ, μ, η)

#maxtime = 368.0
maxtime = 100.0
maxtaxa = 50_000
starting_state = 13
starting_state = 21



trees = Dict()
for seed in 7:100
    Random.seed!(seed)
    tree = sim_bdshift(bds3, maxtime, maxtaxa, starting_state)
    if (ntaxa(tree) > 1)
        println("seed: $seed")
        prune_extinct!(tree)
        begin
            println("ntaxa: $(length(tree.Leaves))")
            Ntotal = sum(sum(branch.N) for (key, branch) in tree.Branches)
            println("N: $Ntotal")
        end
        trees[seed] = tree
    end
end




Random.seed!(68)
tree = sim_bdshift(bds3, maxtime, maxtaxa, starting_state)
prune_extinct!(tree)


nw = newick(tree, bds3)

nw
sizeof(nw)















