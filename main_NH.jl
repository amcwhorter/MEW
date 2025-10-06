include("beano2.2.jl")
using JSON, Graphs, Shapefile, DataFrames, GraphPlot, ProgressBars, Serialization


function main(alph, bet, m1, m2, initialization=nothing)

    data = JSON.parsefile("NH/NH_dual_graph_stripped.json")
    nodes = data["nodes"]
    links = data["links"]

    # is_boundary = [node["boundary_node"] for node in nodes]
    # boundary_lengths = zeros(320)
    # for i in 1:320
        # if is_boundary[i]
            # boundary_lengths[i] = nodes[i]["boundary_perim"]
        # end
    # end

    # perim_dict = Dict()
    g = Graphs.SimpleGraph(length(nodes))
    for link in links
        u, v = link["source"], link["target"]
        e = simple_edge(u + 1, v + 1)
        # perim_dict[e] = link["shared_perim"]
        add_edge!(g, e)
    end

    # areas = [node["area"] for node in nodes]

    table = Shapefile.Table("NH/NH.shp")
    df = DataFrame(table)

    epsilon = 0.1
    if isnothing(initialization) 
        marked_edges = Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}()
        trees = []
        while true
            tree = wilsons(g)
            marked_edge = find_balanced_edge_cut(tree, df)
            if !isnothing(marked_edge)
                push!(marked_edges, marked_edge)
                push!(trees, tree)
                break
            end
        end
        tree = trees[1]
    else
        tree = initialization[1]
        marked_edges = initialization[2]
    end

    # c_x = [node["C_X"] for node in grafton_nodes]
    # c_y = [node["C_Y"] for node in grafton_nodes]
    # p = gplot(tree, c_x, c_y)
    # display(p)

    # p2 = gplot(g, c_x, c_y)
    # display(p2)

    num_iterations = 1000
    accept = 0
    reject = 0
    cycle_inters = 0
    partitions = []
    for i in 1:num_iterations

        tree_old = copy(tree)
        marked_edges_old = copy(marked_edges)
        
        cycle_edges, edge_plus, tree_new, old_edge, new_edge, marked_edges_new, tries = proposal(g, df, tree, marked_edges, epsilon, 2)
        reject += tries

        a1 = transition_probability(cycle_edges, edge_plus, old_edge, new_edge, marked_edges_old, marked_edges_new, tree_old, tree_new)
        a2 = calculate_τs(g, tree_old, tree_new, marked_edges_old, marked_edges_new)

        # a3 = energy_function(g, tree_new, marked_edges_new) / energy_function(g, tree_old, marked_edges_old)
        a3 = multivariate_energy_function(df, g, tree_new, marked_edges_new, alph, bet, m1, m2) / multivariate_energy_function(df, g, tree_old, marked_edges_old, alph, bet, m1, m2)



        a = log(a1) + a2 + log(a3)

        

        if log(rand(Float64)) < a
            accept += 1
            if a1 != 1
                cycle_inters += 1
            end
        else
            tree_new = copy(tree_old)
            marked_edges_new = copy(marked_edges_old)
        end
        push!(partitions, partition(tree_new, marked_edges_new))
        tree = copy(tree_new)
        marked_edges = copy(marked_edges_new)
    end

    return partitions, tree, marked_edges
end

