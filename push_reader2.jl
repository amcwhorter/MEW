using Serialization, Shapefile, DataFrames, Distributions, JSON, ProgressBars
using StatsBase, Graphs, SpecialFunctions
using Distributed, SharedArrays
addprocs(80)

# Make sure all processes have access to the graph and necessary functions
@everywhere begin
    using Serialization, Graphs, StatsBase, DataFrames, Shapefile, JSON, ProgressBars, Distributions, SpecialFunctions
    
    # You'll need to include beano2.2.jl on all processes
    include("beano2.3.jl")
    
    # Function to process a single batch of runs in parallel
    function process_batch_parallel(run_dir, batch_start, batch_end, g, df)
        batch_cuts = []
        batch_mm = []
        
        for i in batch_start:batch_end
            try
                c = deserialize("$(run_dir)/run$(i).jls")
                for j in c[1]
                    push!(batch_cuts, length(cut_edges(j, g)))
                    push!(batch_mm, top_percent_from_partitions(j, df))
                end
            catch e
                println("Error processing $(run_dir)/run$(i).jls: $e")
            end
        end
        
        return batch_cuts, batch_mm
    end
end

# Load your data (keep this part the same)
table = Shapefile.Table("NH/NH.shp")
df = DataFrame(table)
data = JSON.parsefile("NH/NH_dual_graph_stripped.json")
nodes = data["nodes"]
links = data["links"]
g = Graphs.SimpleGraph(length(nodes))
for link in links
    u, v = link["source"], link["target"]
    e = simple_edge(u + 1, v + 1)
    add_edge!(g, e)
end

# Modified process_runs_batch function with pmap parallelization
function process_runs_batch_parallel(run_dirs, batch_size=10, num_runs=100)
    original_shape = size(run_dirs)
    
    all_cuts = Array{Vector{Int}}(undef, original_shape...)
    all_mm = Array{Vector{Float64}}(undef, original_shape...)  # ADDED: Missing declaration
    
    for idx in CartesianIndices(run_dirs)
        run_dir = run_dirs[idx]
        println("Processing $run_dir")
        
        # Create batch ranges
        batch_ranges = []
        for batch_start in 1:batch_size:num_runs
            batch_end = min(batch_start + batch_size - 1, num_runs)
            push!(batch_ranges, (batch_start, batch_end))
        end
        
        # Use pmap to process batches in parallel
        batch_results = pmap(batch_ranges) do (batch_start, batch_end)  # FIXED: Removed extra parenthesis, assigned to variable
            process_batch_parallel(run_dir, batch_start, batch_end, g, df)
        end
        
        # FIXED: Extract cuts and mm from batch_results
        batch_cuts = [result[1] for result in batch_results]
        batch_mm = [result[2] for result in batch_results]
        
        # Combine results
        all_cuts[idx] = vcat(batch_cuts...)
        all_mm[idx] = vcat(batch_mm...)
        
        GC.gc()
    end
    
    return all_cuts, all_mm
end

# Keep your other functions unchanged
function load_cuts_incrementally(run_dir, num_runs)
    cuts = []
    for i in 1:num_runs
        c = deserialize("$(run_dir)/run$(i).jls")
        append!(cuts, c[1])
    end
    return cuts
end

function compute_pairwise_ks(all_cuts)
    n = length(all_cuts)
    ks_matrix = zeros(n, n)
    
    for i in 1:n
        GC.gc()
        data_i = all_cuts[i]
        for j in i+1:n
            data_j = all_cuts[j]
            ks_matrix[i,j] = ks_two_sample(data_i, data_j)
            ks_matrix[j,i] = ks_matrix[i,j]
        end
    end
    
    return ks_matrix
end

# ADDED: Missing build_count_grid function
function build_count_grid(data, xs, ys)
    nx = length(xs)
    ny = length(ys)
    grid = zeros(Int, nx, ny)
    
    for row in eachrow(data)
        x_idx = searchsortedlast(xs, row[1])
        y_idx = searchsortedlast(ys, row[2])
        if x_idx > 0 && y_idx > 0
            for i in x_idx:nx
                for j in y_idx:ny
                    grid[i, j] += 1
                end
            end
        end
    end
    
    return grid
end

function ks_two_sample_2d_sorted(data1, data2)
    # Ensure data is N×2 matrix
    function ensure_matrix(data)
        if !(data isa AbstractMatrix)
            return hcat([[x,y] for (x,y) in data]...)'
        else
            return data
        end
    end

    dat1 = ensure_matrix(data1)
    dat2 = ensure_matrix(data2)
    all_data = vcat(dat1, dat2)
    xs = sort(unique(all_data[:,1]))
    ys = sort(unique(all_data[:,2]))

    function index_in(sorted, val)
        return searchsortedlast(sorted, val)
    end

    grid1 = build_count_grid(dat1, xs, ys)
    grid2 = build_count_grid(dat2, xs, ys)

    n1 = size(dat1,1)
    n2 = size(dat2,1)

    # Evaluate ECDFs at every unique rectangle corner (xs, ys grid)
    maxdist = 0.0
    for i in 1:length(xs)
        for j in 1:length(ys)
            v1 = grid1[i, j]/n1
            v2 = grid2[i, j]/n2
            d = abs(v1 - v2)
            if d > maxdist
                maxdist = d
            end
        end
    end
    return maxdist
end

function empirical_cdf(data)
    sorted_data = sort(data)
    n = length(sorted_data)
    cumulative_probs = collect(1:n) ./ n
    return function(x)
        idx = searchsortedlast(sorted_data, x)
        return idx == 0 ? 0.0 : cumulative_probs[idx]
    end
end

function ks_two_sample(data1, data2)
    ecdf1 = empirical_cdf(data1)
    ecdf2 = empirical_cdf(data2)
    d = 0.0
    for x in data1
        d = max(d, abs(ecdf1(x) - ecdf2(x)))
    end
    return d
end


@everywhere begin
    function compute_single_ks_2_samp_parallel(i, j)
        GC.gc()
        
        current_dirs = ["push_runs2/10_0.1_50_72/run$(k)" for k in 1:10]
        try
            
            cuts, mm = process_runs_batch_parallel(current_dirs, 10, 10000)
            return cuts, mm

        catch e
            println("Error on i = $i, j = $j: $e")
            return NaN
        end
    end
end


cuts, mm = compute_single_ks_2_samp_parallel(1, 1)
using DelimitedFiles
writedlm("cuts.csv", cuts, ',')
writedlm("mm.csv", mm, ',')