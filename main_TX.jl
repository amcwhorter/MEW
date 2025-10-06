using JSON, Graphs, Shapefile, DataFrames, GraphPlot, ProgressBars, Serialization, Plots

include("beano2.3.jl")

function main(alph, bet, m1, m2, initialization=nothing)
    data = JSON.parsefile("Texas/TX_dual_graph_stripped.json")
    nodes = data["nodes"]
    links = data["links"]

    g = Graphs.SimpleGraph(length(nodes))

    for link in links
        u, v = link["source"], link["target"]
        e = simple_edge(u + 1, v + 1)
        add_edge!(g, e)
    end

    table = Shapefile.Table("Texas/tx/tx.shp")
    df = DataFrame(table)
    
    epsilon = 0.05
    if isnothing(initialization)

        marked_edgess = Vector{Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}}()
        trees = []
        while true
            districts = find_k_partition(g, df, 38, "TOTPOP", 0.01)
            if !isnothing(districts)
                tree, marked_edges = partition_to_tree_marked_edges(g, districts)
                push!(trees, tree)
                push!(marked_edgess, marked_edges)
                break
            end
        end
        tree = trees[1]
        marked_edges = marked_edgess[1]
    else 
        tree = initialization[1]
        marked_edges = initialization[2]
    end

    num_iterations = 1000
    accept = 0
    cycle_inters = 0
    partitions = []
    
    a2_tracker = []
    a3_tracker = []

    for i in 1:num_iterations
        tree_old, marked_edges_old = deepcopy(tree), deepcopy(marked_edges)
        cycle_edges, edge_plus, tree_new, old_edge, new_edge, marked_edges_new, tries = proposal(g, df, tree_old, marked_edges_old, 0.1, 38)
        

        a1 = transition_probability(cycle_edges, edge_plus, old_edge, new_edge, marked_edges_old, marked_edges_new, tree_old, tree_new)
        a2 = calculate_τs(g, tree_old, tree_new, marked_edges_old, marked_edges_new)
        a3 = TX_energy_function(df, g, tree_new, marked_edges_new, 0.2, bet, 38, m2) - TX_energy_function(df, g, tree_old, marked_edges_old, 0.2, bet, 38, m2)
        # a3 = 0
        
        a = log(a1) + a2 + a3

        push!(a3_tracker, a3)
        push!(a2_tracker, a2)

        if log(rand(Float64)) < a
            accept += 1
            if a1 !=1
                cycle_inters += 1
            end
        else
            tree_new = copy(tree_old)
            marked_edges_new = deepcopy(marked_edges_old)
        end

        push!(partitions, partition(tree_new, marked_edges_new))
        tree = deepcopy(tree_new)
        marked_edges = deepcopy(marked_edges_new)
    end
    return partitions, tree, marked_edges, a2_tracker, a3_tracker
end

    