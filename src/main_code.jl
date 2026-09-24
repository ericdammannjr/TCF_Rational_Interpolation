using CairoMakie 
include("funs/funs.jl")

# These are the nodes (and corresponding weights) computed by the CGA

tolerance = 1e-12
nodes, errors_history = cga(-1,1; tol = tolerance)
w = eval_weights(nodes, eval_func.(nodes))

# Then add to the computed nodes what would be the test nodes in with more iteration of the CGA (except with 10 nodes between intervals)

test_nodes = Vector{Float64}()

for m in 1:length(nodes)-1

    l = nodes[m]
    r = nodes[m+1]
    h = (r-l)/(10+1)
    append!(test_nodes,[l+k*h for k in 1:10])

end

append!(test_nodes, nodes)
sort!(test_nodes)
rat_approx = eval_tcf.(test_nodes, Ref(nodes), Ref(w))

# This figure shows the absolute error of the rational approximation

fig1 = Figure()

ax1 = Axis(fig1[1,1],
    xlabel = "x",
    ylabel = "y"
)

lines!(ax1, test_nodes, rat_approx .- eval_func.(test_nodes),
    color= :blue,
)

scatter!(ax1, nodes, zeros(length(nodes)),
    markersize = 10,
    color = :red
)

fig1

save("figs/error_plot.png", fig1)

# This figure shows the convergence of the rational approximation

fig2 = Figure()

ax2 = Axis(fig2[1,1],
    xlabel = "Number of nodes",
    ylabel = "Maximum error",
    yscale = log10
)

num_nodes = 3 .+ collect(1:length(errors_history))

lines!(
    ax2,
    num_nodes,
    errors_history,
    color = :blue,
    linewidth = 2
)

scatter!(
    ax2,
    num_nodes,
    errors_history,
    color = :blue,
    markersize = 8
)

fig2

save("figs/convergence_plot.png", fig2)