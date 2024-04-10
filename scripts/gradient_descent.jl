
using Revise
using Pesto
using CairoMakie
using ForwardDiff
using Optim
using ProgressMeter
using Distributions

## todo :
##
## * test on fish tree                          
## * use greedy approach                        
## * figure out error handling/why crash?    
##      * happens in solving D(t)
## * do some time benchmarks                    
## * also do CBDP on Newton method?          

phy = readtree(Pesto.path("primates.tre"))
ρ = 0.635
primates = SSEdata(phy, ρ)

res, model, i = optimize_newton(primates)
logL_root(model, primates)



results = []
for i in 1:40
    res, model, i = optimize_newton(primates)
    push!(results, res)
end



upper = [0.4, 2.0, 1.0]
g, h = logistic(upper, 0.5)


fig = Figure(size = (600, 400))
ax1 = Axis(fig[1,1], xlabel = "λ", ylabel = "μ", xticklabelrotation=45.0)
ax2 = Axis(fig[2,1], xlabel = "λ", ylabel = "η", xticklabelrotation=45.0)
ax3 = Axis(fig[1,2], xlabel = "η", ylabel = "μ", xticklabelrotation=45.0)
ax4 = Axis(fig[1,3], xlabel = "r = λ - μ", ylabel = "logl", xticklabelrotation=45.0, title = "likelihoods")
ax5 = Axis(fig[2,3], xlabel = "replicate", ylabel = "logL", title = "likelihoods")

for (i,res) in enumerate(results)
    x = g(res.initial_x)
    μ = x[1] + x[2]
    λ = x[1] + x[2] + x[3]
    r = λ - μ
    η = x[1]
    logl = -res.minimum

    scatter!(ax1, λ, μ, color = :black)
    scatter!(ax2, λ, η, color = :black)
    scatter!(ax3, η, μ, color = :black)
    scatter!(ax4, r, logl, color = :black)
    scatter!(ax5, i, logl, color = :black)
end

title = Label(fig[0,1:2], "starting values", fontsize = 15)

mlindex = argmax([-res.minimum for res in results])

x = g(results[mlindex].minimizer)
λ = x[1] + x[2] + x[3]
μ = x[1] + x[2]
r = λ - μ
η = x[1]
scatter!(ax1, λ, μ, color = :red)
scatter!(ax2, λ, η, color = :red)
scatter!(ax3, η, μ, color = :red)

#ylims!(ax4, -685, -684)
colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)
colsize!(fig.layout, 1, Relative(0.33))
fig


model = Pesto.newmodel(x)
rates = birth_death_shift(model, primates)

#rates[1:end-1,:mean_netdiv]
hist(rates[1:end-1,:mean_netdiv])
rml

#rates[!,:nshift] 

filter(:shift_bf => x -> x > 3.0, rates)



## try the fish tree

phy = readtree("/home/bkopper/projects/pesto_empirical/data/empirical/Actinopterygii_Rabosky2018.tree")
ρ = 0.37
fishtree = SSEdata(phy, ρ)

results = []
models = []
is = []

for _ in 1:3
    res, model, i = optimize_newton(fishtree)
    push!(results, res)
    push!(models, model)
    push!(is, i)
end


xinit = [0.009496309224430554, 0.0487757521015954, 0.018318692958861803]

res, model, i = optimize_hyperparameters(fishtree; n = 10)

@elapsed logL_root(model, fishtree)





save("/tmp/newton_method_starting_values.pdf", fig)




upper = [0.40, 2.0, 1.0]
g, h = logistic(upper, 0.5)

rml, μml = estimate_constant_netdiv_mu(primates)

dη = Distributions.LogNormal(log(0.01), 0.5)
dμ = Distributions.LogNormal(log(μml), 0.5)
dr = Distributions.LogNormal(log(rml), 0.5)

## truncate the distribution
dη = Distributions.Truncated(dη, 0.0, upper[1])
dμ = Distributions.Truncated(dμ, 0.0, upper[2])
dr = Distributions.Truncated(dr, 0.0, upper[3])

for i in 1:100_000
    xinit = zeros(3)
    xinit[1] = rand(dη)
    xinit[2] = rand(dμ)
    xinit[3] = rand(dr)

    println(xinit)
    xinit_tilde = h(xinit)
end



d = Distributions.LogNormal(log(1.0), 0.5)
d_truncated = Distributions.Truncated(d, 0.0, 1.0)

sum(rand(d_truncated, 10000) .> 1.0)











x0 = [rand(d) for d in (d1, d2, d3)]
x0 .< upper
x0_tilde = h(x0)
res = Optim.optimize(f, g!, h!, x0_tilde, inner_optimizer)

model = Pesto.newmodel(x0)

logL_root(model, primates)



x = collect(range(-10, 10; length = 100))

plot(x, g.(x))

g(0.8)









optres, model, i, history = optimize_hyperparameters(primates, n_attempts = 6)

results = optimize_hyperparameters(primates, n_attempts = 20, n_start_positions = 5)



λ = [x[1] for x in history]
μ = [x[2] for x in history]
η = [x[3] for x in history]


fig = Figure(size = (800, 800))
ax1 = Axis(fig[1,1], xlabel = "iteration", ylabel = "λ", yscale = log10)
lines!(ax1, eachindex(history), λ)
ax2 = Axis(fig[1,2], xlabel = "iteration", ylabel = "μ", yscale = log10)
lines!(ax2, eachindex(history), μ)
ax3 = Axis(fig[2,1], xlabel = "iteration", ylabel = "η", yscale = log10)
lines!(ax3, eachindex(history), η)
ax4 = Axis(fig[2,2], xlabel = "μ", ylabel = "λ", yscale = log10, xscale = log10)
scatter!(ax4, μ, λ, color = (:black, 0.05))
ax5 = Axis(fig[3,1], xlabel = "λ", ylabel = "η", yscale = log10, xscale = log10)
scatter!(ax5, λ, η, color = (:black, 0.05))
ax6 = Axis(fig[3,2], xlabel = "μ", ylabel = "η", yscale = log10, xscale = log10)
scatter!(ax6, μ, η, color = (:black, 0.05))
fig

models = [Pesto.newmodel(result.minimizer) for result in results]

ratesx = [birth_death_shift(model, primates) for model in models]
logls = [logL_root(model, primates) for model in models]


treeplot(primates, ratesx[3])
treeplot(primates, rates2)

rates1[!,:nshift] |> sum
rates2[!,:nshift] |> sum

argmax(logls)

[model.η for model in models]

[x.f_calls for x in results]

x1 = results[1]

x1.
















