using Glob
using CairoMakie
using CSV, DataFrames
using Statistics


## load simulated data
fpaths = Glob.glob("output/simulations/age_scaling_effect/shift_rate_through_time/*.csv")
dfs_simulated = []
for fpath in fpaths
    push!(dfs_simulated, CSV.read(fpath, DataFrame))
end

## load empirical data
names_selected = [
    "Actinopterygii", "Mammalia", "Rosidae", "Chondrichthyes", "Squamata", "Asteraceae", 
    "Agaricomycetes", "Lissamphibia", "Aves", "Polypodiophyta", "Lecanoromycetes", "Cichlidae"
    ]
fpaths = Glob.glob("output/empirical_joint/shift_rate_through_time/*.csv")
dfs_empirical = []
dfs_empirical_all = []
names = []
names_all = []
for fpath in fpaths
    name = split(Base.basename(fpath), "_")[1]
    push!(dfs_empirical_all, CSV.read(fpath, DataFrame))
    push!(names_all, name)
    
    if name in names_selected
        if name == "Polypodiophyta"
            if !occursin("Lehtonen", fpath)
                continue
            end
        end
        if name == "Agaricomycetes"
            if !occursin("SanchezGarcia", fpath)
                continue
            end
        end
        push!(dfs_empirical, CSV.read(fpath, DataFrame))
        push!(names, name)
    end
end

dfs_empirical
names








fig = Figure(size = (750, 650));


ax_summary = Axis(fig[1:3,1:2],
                xgridvisible = false, 
                ygridvisible = false,
                topspinevisible = false,
                rightspinevisible = false,
                yticks = ([1.0, 2.0], [L"\text{increasing}", L"\text{decreasing}"]),
                yticklabelrotation=π/2,
                xlabel = L"\text{frequency}",
                xticklabelsize = 9)
                #yticklabelsize = 9)
Label(fig[0,1:2],L"\text{a) summary}")


axs_simulation = []
Label(fig[0,4:7],L"\text{b) simulated phylogenies}")

for i in 1:3
    for j in 1:3
        ax = Axis(fig[j, i+3], xreversed = true,
                xgridvisible = false, 
                ygridvisible = false,
                topspinevisible = false,
                rightspinevisible = false,
                xticklabelsize = 9,
                yticklabelsize = 9)
        push!(axs_simulation, ax)
    end
end

axs_empirical = []
Label(fig[4,1:6],L"\text{c) empirical phylogenies}")

q = [1]
for j in 1:6
    for i in 1:2
        ax = Axis(fig[4+i, j], xreversed = true,
                title = names[q[1]],
                xgridvisible = false, 
                ygridvisible = false,
                topspinevisible = false,
                rightspinevisible = false,
                titlesize = 10,
                xticklabelsize = 9,
                yticklabelsize = 9)
        q[1] += 1

        push!(axs_empirical, ax)
    end
end


## add lines
for (ax, df) in zip(axs_simulation, dfs_simulated)
    lines!(ax, df[!,:times], df[!,:nshift_fitted], color = :gray)
    lines!(ax, df[!,:times], df[!,:nshift_true], color = :black)
end

for (ax, df) in zip(axs_empirical, dfs_empirical)
    lines!(ax, df[!,:times], df[!,:nshift], color = :orange)
end

## add x-axis label
Label(fig[7,1:6], L"\text{time before present (Ma)}")

## add summary barplot
present_minus_mrca_sim_true = [df[!,:nshift_true][end] - df[!,:nshift_true][1] for df in dfs_simulated]
present_minus_mrca_sim_fitted = [df[!,:nshift_fitted][end] - df[!,:nshift_fitted][1] for df in dfs_simulated]
present_minus_mrca_empirical = [df[!,:nshift][end] - df[!,:nshift][1] for df in dfs_empirical_all]


freq_increasing_sim_true = mean(present_minus_mrca_sim_true .> 0) 
freq_increasing_sim_fitted = mean(present_minus_mrca_sim_fitted .> 0) 
freq_increasing_empirical = mean(present_minus_mrca_empirical .> 0) 

tbl = (
    height = [freq_increasing_sim_true, 1.0 - freq_increasing_sim_true, 
              freq_increasing_sim_fitted, 1.0 - freq_increasing_sim_fitted,
              freq_increasing_empirical, 1.0 - freq_increasing_empirical ],
    cat = [1,2,1,2,1,2],
    grp = [1,1,2,2,3,3]
)

colors = [:black, :gray, :orange]

barplot!(ax_summary, 
    tbl.cat,
    tbl.height,
    dodge = tbl.grp,
    color = colors[tbl.grp],
    direction = :x
)

# Legend
labels = ["simulated (true par.)", "simulated (estimated par.) ", "empirical (estimated par.)"]
elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
title = ""
Legend(fig[1,2], elements, labels, title, tellheight = true, tellwidth = true, patchsize = (10, 10), labelsize = 7,
        framevisible = false)


for i in 1:6
    colsize!(fig.layout, i, Relative(1/6.5))
end
colgap!(fig.layout, 4)

for i in [1,2,3,5,6]
    rowsize!(fig.layout, i, Relative(1/5.5))
end
rowgap!(fig.layout, 6)

#Label(fig[1:3,0], L"\text{frequency}", rotation = π/2)
Label(fig[5:6,0], L"\text{mean shift rate }(\frac{dN}{dt}(t))", rotation = π/2)
Label(fig[1:3,3], L"\text{mean shift rate }(\frac{dN}{dt}(t))", rotation = π/2)

for i in 1:6
    colsize!(fig.layout, i, Relative(1/6.33))
end
colsize!(fig.layout,0, Relative(1/18))

fig

CairoMakie.save("figures/shift_through_time.pdf", fig)