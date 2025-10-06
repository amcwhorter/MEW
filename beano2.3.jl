using Graphs, GraphPlot, StatsBase, Combinatorics, LinearAlgebra, SparseArrays
#contains helpful functions - log metropolis hastings version

function simple_edge(edge::Graphs.SimpleGraphs.SimpleEdge)
    u, v = src(edge), dst(edge)
    if u < v
        return edge
    else 
        return Graphs.SimpleGraphs.SimpleEdge(v, u)
    end
end

function simple_edge(u::Int, v::Int)
    if u < v
        return Graphs.SimpleGraphs.SimpleEdge(u, v)
    else
        return Graphs.SimpleGraphs.SimpleEdge(v, u)
    end
end



function instantiate(m::Int, n::Int)
    g = Graphs.grid((m, n))
    min_tree_edges = kruskal_mst(g, Graphs.weights(g))
    tree = SimpleGraph(m * n)

    for edge in min_tree_edges 
        add_edge!(tree, edge)
    end

    tree = SimpleGraph(tree)
    marked_edges = sample(min_tree_edges, 1, replace=false)

    return g, tree, marked_edges
end


function instantiate(g::Graphs.SimpleGraph)
    min_tree_edges = kruskal_mst(g, Graphs.weights(g))
    tree = SimpleGraph(length(collect(1:nv(g))))

    for edge in min_tree_edges
        add_edge!(tree, simple_edge(edge))
    end

    tree = SimpleGraph(tree)
    marked_edges = sample(min_tree_edges, 4, replace=false)
    
    return g, tree, marked_edges
end


function grid_layout(g)
    xlocs = Vector{Int64}()
    ylocs = Vector{Int64}()
    for i in 1:m
        for j in 1:n
            push!(xlocs, i)
            push!(ylocs, j)
        end
    end
    return xlocs, ylocs
end

function tally(df, column_name, partition_dict)
    to_add = df[:, column_name]

    totals = zeros(maximum(values(partition_dict)))

    for i in keys(partition_dict)
        j = partition_dict[i]
        totals[j] += to_add[i]
    end 
    return totals
end

function within_percent_of_ideal(vals, ideal, epsilon)
    ideals = ideal * ones(length(vals))
    epsilons = abs.((vals - ideals)/ideal)
    return maximum(epsilons) < epsilon
end

function find_balanced_edge_cut(tree, df)
    for edge in edges(tree)
        part_dict = partition(tree, [edge])

        totals = tally(df, "TOTPOP", part_dict)

        ideal = sum(df[:,"TOTPOP"])/2

        if within_percent_of_ideal(totals, ideal, 0.01)
            return edge
        end
    end

    println("no balanced edge cuts")

    return 

end

function cycle_basis_step(g, tree, marked_edges)
    edges_g = Set(edges(g))
    edges_tree = Set(edges(tree))
    difference_set = setdiff(edges_g, edges_tree)

    edge_plus = rand(difference_set)
    add_edge!(tree, edge_plus)

    cycle = Graphs.cycle_basis(tree)[1]
    cycle_edges = simple_edge.(cycle, circshift(cycle, -1))

    possible_cuts = setdiff(cycle_edges, marked_edges)
    edge_minus = rand(possible_cuts)

    rem_edge!(tree, edge_minus)

    return cycle_edges, edge_plus
end

function marked_edge_step(tree, marked_edges)
    old_edge = rand(marked_edges)
    v1, v2 = src(old_edge), dst(old_edge)
    chosen_v = rand((v1, v2))
    new_edge = simple_edge(chosen_v, rand(Graphs.neighbors(tree, chosen_v)))
    setdiff!(marked_edges, [old_edge])
    push!(marked_edges, new_edge)
    return old_edge, new_edge
end



function transition_probability(cycle_edges, edge_plus, marked_edge_old, marked_edge_new, marked_edges_old, marked_edges_new, tree_old, tree_new)
    if marked_edge_new == edge_plus
        return 0
    end

    w, x = src(marked_edge_old), dst(marked_edge_old)
    u, v = src(marked_edge_new), dst(marked_edge_new)

    if Set([u, v]) == Set([w, x])
        d_u = length(neighbors(tree_old, u))
        d_u_p = length(neighbors(tree_new, u))
        d_v = length(neighbors(tree_old, v))
        d_v_p = length(neighbors(tree_new, v))

        pm = (d_u + d_v)/(d_u_p + d_v_p) * d_u_p/d_u * d_v_p/d_v
    else
        u = intersect([u, v], [w, x])[1]
        d_u = length(neighbors(tree_old, u))
        d_u_p = length(neighbors(tree_new, u))
        pm = d_u_p/d_u
    end


    cycle = Set(cycle_edges)
    m_old = Set(marked_edges_old)
    m_new = Set(marked_edges_new)

    l = length(setdiff(cycle, m_old))
    l_p = length(setdiff(cycle, m_new))

    pt = l/l_p

    return pt * pm
end

function remove_zero_rows_and_cols(lapmat)
    # Identify zero rows where the number of non-zero elements is zero
    zero_rows = findall(row -> count(!iszero, lapmat[row, :]) == 0, 1:size(lapmat, 1))

    # Determine indices for non-zero rows and columns
    non_zero_indices = setdiff(1:size(lapmat, 1), zero_rows)

    # Subset the matrix to exclude zero rows and corresponding columns
    cleaned_lapmat = lapmat[non_zero_indices, non_zero_indices]

    return cleaned_lapmat
end


function log_τ(edgelist)

    g = SimpleGraphFromIterator(edgelist)

    L = laplacian_matrix(g)

    L = Matrix(remove_zero_rows_and_cols(L))
    
    Λ = cholesky(L[2:end,2:end])
    U = Λ.L
    c = 2 * sum(log.([U[i, i] for i in 1:size(U)[1]]))

    return c
end

function nodes_to_edgelist(g, part)
    connected_edges = [simple_edge(u, v) for (u, v) in combinations(part, 2)]
    part_edges = intersect(connected_edges, edges(g))
    return part_edges
end

function cut_edges(ptition, g)
    cuts = []
    for e in edges(g)
        u, v = src(e), dst(e)

        if ptition[u] != ptition[v]
            push!(cuts, simple_edge(u, v))
        end
    end
    return cuts
end

# function cut_edges(g, tree, marked_edges)

#     for m in marked_edges
#         rem_edge!(tree, m)
#     end

#     parts = connected_components(tree)
#     es = [nodes_to_edgelist(g, part) for part in parts]

#     cuts = setdiff(edges(g), Set(Iterators.flatten(es)))

#     for m in marked_edges
#         add_edge!(tree, m)
#     end

#     return length(cuts)
# end

function partition(parts::Vector{Vector{Int}})
    partition_dict = Dict{Int, Int}()
    for (i, part) in enumerate(parts)
        for node in part
            partition_dict[node] = i
        end
    end
    return partition_dict
end 


function partition(tree::SimpleGraph{Int64}, marked_edges::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}})
    for m in marked_edges
        rem_edge!(tree, m)
    end

    parts = connected_components(tree)

    for m in marked_edges
        add_edge!(tree, m)
    end

    return partition(parts)
end



function borders_dict(cut_edges, p_dict)
    borders_dict = Dict()
    for edge in cut_edges
        u, v = src(edge), dst(edge)
        key = Set([p_dict[u], p_dict[v]])
        if !(key in keys(borders_dict))
            borders_dict[key] = 1
        else
            borders_dict[key] += 1
        end
    end
    return borders_dict
end


function calculate_τs(g, tree_old, tree_new, marked_edges_old, marked_edges_new)
    for m in marked_edges_old
        rem_edge!(tree_old, m)
    end

    for m in marked_edges_new
        rem_edge!(tree_new, m)
    end


    #calculate τ relating to border lengths
    old_parts = connected_components(tree_old)
    new_parts = connected_components(tree_new)

    p_dict_old = partition(old_parts)
    p_dict_new = partition(new_parts)

    cut_edges_old = cut_edges(p_dict_old, g)
    cut_edges_new = cut_edges(p_dict_new, g)

    borders_dict_old = borders_dict(cut_edges_old, p_dict_old)
    borders_dict_new = borders_dict(cut_edges_new, p_dict_new)

    τ_border_lengths = sum(log.(values(borders_dict_old))) - sum(log.(values(borders_dict_new)))

    #calculate τ relating to spanning trees
    # for part in copy(old_parts)
    #     for part_ in copy(new_parts)
    #         if Set(part) == Set(part_)
    #             setdiff!(old_parts, [part])
    #             setdiff!(new_parts, [part_])
    #         end
    #     end
    # end

    old_edges = [nodes_to_edgelist(g, part) for part in old_parts]
    new_edges = [nodes_to_edgelist(g, part) for part in new_parts]

    old_τs = []
    for edgelist in old_edges
        push!(old_τs, log_τ(edgelist))
    end

    new_τs = []
    for edgelist in new_edges
        push!(new_τs, log_τ(edgelist))
    end

    # old_τs = [log_τ(edgelist) for edgelist in old_edges]
    # new_τs = [log_τ(edgelist) for edgelist in new_edges]

    τ_spanning_trees = sum(old_τs) - sum(new_τs)
    
    τ = τ_border_lengths + τ_spanning_trees

    for m in marked_edges_old
        add_edge!(tree_old, m)
    end

    for m in marked_edges_new
        add_edge!(tree_new, m)
    end

    return τ
end


function energy_function(g, tree, marked_edges, beta)

    ptition = partition(tree, marked_edges)
    cuts = length(cut_edges(ptition, g))

    J = (12 - cuts) ^ 2 
    
    return exp(-beta * J)
end

function energy_function(a::Number)
    J = (12 - a) ^ 2

    return exp(- beta * J)
end

function wilsons(g)
    spanning_tree = Set()

    visited = Set()

    start = rand(1:nv(g))
    push!(visited, start)

    while length(visited) < nv(g)
        unvisited = setdiff(1:nv(g), visited)
        current = rand(unvisited)
        path = [current]

        while current ∉ visited

            next = rand(all_neighbors(g, current))

            loop_index = findfirst(==(next), path)
            
            if !isnothing(loop_index)
                path = path[1:loop_index]
            else
                push!(path, next)
            end

            current = next
        end

        for i in 1:(length(path) - 1)
            u, v = path[i], path[i + 1]
            push!(visited, u)
            push!(spanning_tree, simple_edge(u, v))
        end
    end

    tree = Graphs.SimpleGraph(nv(g))
    
    for edge in spanning_tree
        add_edge!(tree, edge)
    end


    return tree

end
