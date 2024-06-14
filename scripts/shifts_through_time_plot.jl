using Glob
using CairoMakie


## load simulated data
fpaths = Glob.glob("output/simulations/age_scaling_effect/shift_rate_through_time/*.csv")
dfs_simulated = []
for fpath in fpaths
    push!(dfs_simulated, CSV.read(fpath, DataFrame))
end

## load empirical data
fpaths = Glob.glob("output/empirical_joint/shift_rate_through_time/*.csv")
dfs_empirical = []
names = []
for fpath in fpaths
    push!(dfs_empirical, CSV.read(fpath, DataFrame))
    name = split(Base.basename(fpath), "_")[1]
    push!(names, name)
end



fig = Figure(size = (750, 550));


ax_summary = Axis(fig[1:2,1:2],
                xgridvisible = false, 
                ygridvisible = false,
                topspinevisible = false,
                rightspinevisible = false,
                xticks = ([1.0, 2.0], [L"\text{increasing}", L"\text{decreasing}"]),
                xticklabelsize = 9,
                yticklabelsize = 9)
Label(fig[0,1:2],L"\text{summary}")


axs_simulation = []
Label(fig[0,3:6],L"\text{simulated phylogenies}")

for i in 1:4
    for j in 1:2
        ax = Axis(fig[j, i+2], xreversed = true,
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
Label(fig[3,1:6],L"\text{empirical phylogenies}")

q = [1]
for j in 1:6
    for i in 1:2
        ax = Axis(fig[3+i, j], xreversed = true,
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
Label(fig[6,1:6], L"\text{time before present (Ma)}")

## add summary barplot
present_minus_mrca_sim_true = [df[!,:nshift_true][end] - df[!,:nshift_true][1] for df in dfs_simulated]
present_minus_mrca_sim_fitted = [df[!,:nshift_fitted][end] - df[!,:nshift_fitted][1] for df in dfs_simulated]
present_minus_mrca_empirical = [df[!,:nshift][end] - df[!,:nshift][1] for df in dfs_empirical]

using Statistics

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
    color = colors[tbl.grp]
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

Label(fig[1:2,0], L"\text{frequency}", rotation = π/2)
Label(fig[4:5,0], L"\text{mean shift rate }(\frac{dN}{dt}(t))", rotation = π/2)

for i in 1:6
    colsize!(fig.layout, i, Relative(1/6.33))
end
colsize!(fig.layout,0, Relative(1/18))
