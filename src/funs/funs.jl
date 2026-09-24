function eval_func(z)

    f = log(1+eps(1.0)+z)

    return f

end

# This is a function to compute the TCF weights 

function eval_weights(z_nodes, f)

    n = length(z_nodes)
    w = zeros(eltype(z_nodes), n)
    w[1] = f[1]

    for k = 2:n

        t = f[k]

        for i = 1:k-1

            t = (z_nodes[k] - z_nodes[i]) / (-w[i] + t)

        end

        w[k] = t

    end

    return w

end

# This is a function that computes the TCF for some scalar or vector z

function eval_tcf(z, z_nodes, w)

    n = length(w)
    t = w[end]

    for k = n-1:-1:1

        t = w[k] + (z - z_nodes[k]) / t

    end

    return t

end

# This function computes the nodes for the TCF by a continuum greedy algorithm (CGA)


function cga(a, b; tol = 1e-12, n = 3, maxiter = 200)

    x_test = range(a, b, length = 5000)
    nodes = [a, (a+b)/2, b]
    errors_history = Float64[]

    for i in 1:maxiter
        
        t_nodes = Vector{Float64}()

        for j in 1:length(nodes)-1

            l = nodes[j]
            r = nodes[j+1]
            h = (r-l)/(n+1)
            append!(t_nodes, [l+k*h for k in 1:n])

        end
        
        w = eval_weights(nodes, eval_func.(nodes))

        # Error on the fixed test grid
        
        approx = eval_tcf.(x_test, Ref(nodes), Ref(w))
        true_values = eval_func.(x_test)

        max_error = maximum(abs.(true_values .- approx))
        push!(errors_history, max_error)

        # Use the CGA test points to determine the next node

        cga_errors = abs.(
            eval_func.(t_nodes) .-
            eval_tcf.(t_nodes, Ref(nodes), Ref(w))
        )

        _, index = findmax(cga_errors)

        new_node = t_nodes[index]
        push!(nodes, new_node)
        sort!(nodes)

        if max_error < tol

            println("Converged to tolerance with $(length(nodes)) nodes")
            break

        end

    end

    return nodes, errors_history

end